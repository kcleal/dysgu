# cython: language_level=3

from __future__ import absolute_import
import datetime
import time
import os
import heapq
import logging
import numpy as np
import random
from collections import Counter
import pysam
from sys import stdout
import pandas as pd
from dysgu import coverage, graph, call_component, consensus, io_funcs, re_map, post_call, merge_svs, extra_metrics
from dysgu.map_set_utils cimport EventResult, Py_SimpleGraph
from dysgu.map_set_utils import to_dict, echo
from dysgu import sites_utils
from dysgu.io_funcs import intersecter
import pickle
import gc
import itertools
import multiprocessing
from scipy import stats
from libcpp.vector cimport vector

ctypedef EventResult EventResult_t

np.random.seed(0)
random.seed(0)


def filter_potential(input_events, tree, regions_only):
    potential = []
    cdef EventResult_t i
    for i in input_events:
        if i.site_info:
            potential.append(i)
            continue
        if i.sqc is None:
            i.sqc = -1
        if i.svtype == "INS" and i.svlen_precise == 0 and not i.contig and not i.contig2:
            continue
        # Remove events for which both ends are in --regions but no contig was found
        posA_intersects = intersecter(tree, i.chrA, i.posA, i.posA + 1)
        posB_intersects = intersecter(tree, i.chrB, i.posB, i.posB + 1)
        # Remove events for which neither end is in --regions (if --regions provided)
        if tree and regions_only:
            if not posA_intersects and not posB_intersects:
                continue
        potential.append(i)
    return potential


def component_job(infile, component, regions, event_id, clip_length, insert_med, insert_stdev, insert_ppf, min_supp, lower_bound_support,
                  merge_dist, regions_only, assemble_contigs, rel_diffs, diffs, min_size, max_single_size,
                  sites_index, paired_end, length_extend, divergence, hp_tag):
    potential_events = []
    grp_id = event_id
    cdef EventResult_t event
    for event in call_component.call_from_block_model(infile,
                                                      component,
                                                      clip_length,
                                                      insert_med,
                                                      insert_stdev,
                                                      insert_ppf,
                                                      min_supp,
                                                      lower_bound_support,
                                                      assemble_contigs,
                                                      max_single_size,
                                                      sites_index,
                                                      paired_end,
                                                      length_extend,
                                                      divergence,
                                                      hp_tag):
        if event:
            event.grp_id = grp_id
            event.event_id = event_id
            if event.chrA is not None:
                event.chrA = infile.get_reference_name(event.chrA)
            if event.chrB is not None:
                event.chrB = infile.get_reference_name(event.chrB)
            potential_events.append(event)
            event_id += 1
    if not potential_events:
        return potential_events, event_id
    potential_events = filter_potential(potential_events, regions, regions_only)
    return potential_events, event_id


def process_job(msg_queue, args):
    job_path, infile_path, bam_mode, ref_path, regions_path, clip_length, insert_median, insert_stdev, insert_ppf, min_support, \
    lower_bound_support, merge_dist, regions_only, assemble_contigs, rel_diffs, diffs, min_size,\
        max_single_size, sites_index, paired_end, length_extend, divergence, hp_tag = args
    regions = io_funcs.overlap_regions(regions_path)
    completed_file = open(job_path[:-3] + "done.pkl", "wb")
    pysam.set_verbosity(0)
    infile = pysam.AlignmentFile(infile_path, bam_mode, threads=1, reference_filename=None if bam_mode != "rc" else ref_path)
    pysam.set_verbosity(3)
    event_id = 0
    while 1:
        msg = msg_queue.recv()
        if msg == 0:
            break
        res = msg
        event_id += 1
        potential_events, event_id = component_job(infile, res, regions, event_id, clip_length,
                                                   insert_median,
                                                   insert_stdev,
                                                   insert_ppf,
                                                   min_support,
                                                   lower_bound_support,
                                                   merge_dist,
                                                   regions_only,
                                                   assemble_contigs,
                                                   rel_diffs=rel_diffs, diffs=diffs, min_size=min_size,
                                                   max_single_size=max_single_size, sites_index=sites_index,
                                                   paired_end=paired_end, length_extend=length_extend, divergence=divergence,
                                                   hp_tag=hp_tag)
        if potential_events:
            for res in potential_events:
                pickle.dump(res, completed_file)
    completed_file.close()


def pipe1(args, infile, kind, regions, ibam, ref_genome, sample_name, bam_iter=None):
    procs = args['procs']
    low_mem = args['low_mem']
    tdir = args["working_directory"]
    regions_only = False if args["regions_only"] == "False" else True
    paired_end = int(args["paired"] == "True")
    assemble_contigs = int(args["contigs"] == "True")
    if args["max_cov"] == "auto":
        if args["ibam"] is not None:
            mc = coverage.auto_max_cov(args["max_cov"], args["ibam"])
        else:
            mc = coverage.auto_max_cov(args["max_cov"], args["sv_aligns"])
        args["max_cov"] = mc
    else:
        args["max_cov"] = int(args["max_cov"])
    if not args["ibam"]:
        cov_track_path = tdir
    else:
        cov_track_path = None

    if os.path.exists(os.path.join(tdir, "n_aligned_bases.txt")):
        n_aligned_bases_file = int(open(os.path.join(tdir, "n_aligned_bases.txt"), "r").readline().strip())
        assert n_aligned_bases_file > 0
        find_n_aligned_bases = False
    else:
        find_n_aligned_bases = True

    genome_scanner = coverage.GenomeScanner(infile, args["mq"], args["max_cov"], args["regions"], procs,
                                            args["buffer_size"], regions_only,
                                            kind == "stdin",
                                            clip_length=args["clip_length"],
                                            min_within_size=args["min_size"],
                                            cov_track_path=cov_track_path,
                                            paired_end=paired_end,
                                            bam_iter=bam_iter)
    insert_median, insert_stdev, read_len = -1, -1, -1
    max_dist = None
    if args["template_size"] != "":
        try:
            insert_median, insert_stdev, read_len = list(map(int, args["template_size"].split(",")))
        except:
            raise ValueError("Template-size must be in the format 'INT,INT,INT', for insert-median, insert-stdev, read-length")

    if paired_end:
        insert_median, insert_stdev = genome_scanner.get_read_properties(args["max_tlen"], insert_median, insert_stdev, read_len, ibam)
        if insert_median == -1:
            return [], None
        read_len = genome_scanner.approx_read_length
        max_dist = int(insert_median + (insert_stdev * 5))  # 5
        max_clust_dist = 1 * (int(insert_median + (5 * insert_stdev)))
        if args["merge_dist"] is None:
            args["merge_dist"] = max_clust_dist
        logging.info(f"Max clustering dist: {max_clust_dist}")
        args["divergence"] = 1
    else:
        if args["divergence"] == "auto":
            divergence, _ = genome_scanner.get_read_properties(args["max_tlen"], insert_median, insert_stdev, read_len, ibam, find_divergence=True)
            args["divergence"] = divergence
        else:
            args["divergence"] = float(args["divergence"])
            logging.info(f"Sequence divergence: {args['divergence']}")
        if args["mode"] == "pacbio-sequel2":
            max_dist, max_clust_dist = 35, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 2000
        elif args["mode"] == "pacbio-revio":
            max_dist, max_clust_dist = 50, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 2000
        elif args["mode"] == "nanopore-r9" or args["mode"] == "nanopore-r10":
            max_dist, max_clust_dist = 100, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 2000

    # set upper bound on single-partition size
    max_single_size = min(max(args["max_cov"] * 50, 10000), 100000)  # limited between 5000 - 50,000 reads
    event_id = 0
    block_edge_events = []

    clip_length = args["clip_length"]
    merge_dist = args["merge_dist"]
    min_size = args["min_size"]
    length_extend = args["length_extend"]
    divergence = args["divergence"]
    max_divergence = divergence * 13

    read_buffer = genome_scanner.read_buffer
    sites_info = sites_utils.vcf_reader(args["sites"], infile, args["parse_probs"], sample_name, args["ignore_sample_sites"] == "True", args["sites_prob"], args["sites_pass_only"] == "True")

    cdef Py_SimpleGraph G
    G, node_to_name, bad_clip_counter, sites_adder, n_aligned_bases, hp_tag_found = graph.construct_graph(genome_scanner,
                                            infile,
                                            max_dist=max_dist,
                                            clustering_dist=max_clust_dist,
                                            minimizer_dist=150,
                                            minimizer_support_thresh=args["z_depth"],
                                            minimizer_breadth=args["z_breadth"],
                                            k=12,
                                            m=6,
                                            clip_l=clip_length,
                                            min_sv_size=min_size,
                                            procs=1,
                                            mapq_thresh=args["mq"],
                                            debug=None,
                                            paired_end=paired_end,
                                            read_length=read_len,
                                            contigs=args["contigs"],
                                            norm_thresh=args["dist_norm"],
                                            spd_thresh=args["spd"],
                                            mm_only=args["regions_mm_only"] == "True",
                                            sites=sites_info,
                                            trust_ins_len=args["trust_ins_len"] == "True",
                                            low_mem=low_mem,
                                            temp_dir=tdir,
                                            find_n_aligned_bases=find_n_aligned_bases,
                                            position_distance_thresh=args['sd'],
                                            max_search_depth=args['search_depth'],
                                            max_divergence=max_divergence,
                                            no_phase=args['no_phase']
    )

    sites_index = None
    if sites_adder:
        sites_index = sites_adder.sites_index
    logging.info("Graph constructed")

    if not args['no_phase'] and hp_tag_found:
        logging.info("Using HP tag")

    auto_support = False
    if args["min_support"] != "auto":
        args["min_support"] = int(args["min_support"])
        min_support = args["min_support"]
        logging.info(f"Minimum support: {args['min_support']}")
    else:
        auto_support = True
        genome_length = sum(infile.lengths)
        if find_n_aligned_bases:
            bases = n_aligned_bases
        else:
            bases = n_aligned_bases_file
        assert bases > 0
        cov_estimate = bases / genome_length
        min_support = round(1.5 + 0.05 * cov_estimate)
        args["min_support"] = min_support
        logging.info(f"Inferred minimum support: {min_support}")

    if args["pl"] == "pe":  # reads with internal SVs can be detected at lower support
        lower_bound_support = min_support - 1 if min_support - 1 > 1 else 1
    else:
        lower_bound_support = 2

    component_path = f"{tdir}/components.bin"
    cdef bytes cmp_file = component_path.encode("ascii")  # write components to file if low-mem used
    cdef vector[int] cmp
    G.connectedComponents(cmp_file, low_mem, cmp)
    cdef int length_components = cmp.size()
    if low_mem:
        cmp_mmap = np.memmap(component_path, dtype=np.int32, mode='r')
        length_components = len(cmp_mmap)
    if insert_median != -1:
        insert_ppf = stats.norm.ppf(0.05, loc=insert_median, scale=insert_stdev)
        if insert_ppf < 0:
            insert_ppf = 0
    else:
        insert_ppf = -1
    if paired_end:
        rel_diffs = False
        diffs = 15
    else:
        rel_diffs = True
        diffs = 0.15
    write_index = None
    minhq = None
    consumers = []
    msg_queues = []
    if procs > 1:
        write_index = itertools.cycle(range(procs))
        minhq = [(0, i) for i in range(procs)]
        msg_queues = [multiprocessing.Pipe(duplex=False) for _ in range(procs)]
        for n in range(procs):
            job_path = f"{tdir}/job_{n}.pkl"
            proc_args = ( job_path, args["sv_aligns"], args["bam_mode"], args["reference"], args["regions"], clip_length, insert_median, insert_stdev,
                insert_ppf, min_support, lower_bound_support, merge_dist, regions_only, assemble_contigs,
                rel_diffs, diffs, min_size, max_single_size, sites_index, paired_end, length_extend, divergence, hp_tag_found )
            p = multiprocessing.Process(target=process_job, args=(msg_queues[n][0], proc_args,), daemon=True)
            p.start()
            consumers.append(p)

    num_jobs = 0
    completed = 0
    components_seen = 0
    if procs == 1 and low_mem:
        completed_file = open(f"{tdir}/job_0.done.pkl", "wb")
    else:
        completed_file = None
    cdef int last_i = 0
    cdef int start_i, end_i, ci, cmp_idx, item_index, item  # cmp is a flat array of indexes. item == -1 signifies end of component
    for item_idx in range(length_components):
        if low_mem:
            item = cmp_mmap[item_idx]
        else:
            item = cmp[item_idx]
        if item == -1:
            components_seen += 1
            start_i = last_i
            end_i = item_idx
            last_i = item_idx + 1
            component = np.zeros(end_i - start_i)
            ci = 0
            for cmp_idx in range(start_i, end_i):
                if low_mem:
                    component[ci] = cmp_mmap[cmp_idx]
                else:
                    component[ci] = cmp[cmp_idx]
                ci += 1
            if len(component) > 100_000:
                reduced = graph.break_large_component(G, component, min_support)
                for cmp2 in reduced:
                    res = graph.proc_component(node_to_name, cmp2, read_buffer, infile, G, lower_bound_support,
                                               procs, paired_end, sites_index)
                    if not res:
                        continue
                    event_id += 1
                    if procs == 1:
                        potential_events, event_id = component_job(infile, res, regions, event_id, clip_length,
                                                                   insert_median,
                                                                   insert_stdev,
                                                                   insert_ppf,
                                                                   min_support,
                                                                   lower_bound_support,
                                                                   merge_dist,
                                                                   regions_only,
                                                                   assemble_contigs,
                                                                   rel_diffs=rel_diffs, diffs=diffs,
                                                                   min_size=min_size,
                                                                   max_single_size=max_single_size,
                                                                   sites_index=sites_index,
                                                                   paired_end=paired_end,
                                                                   length_extend=length_extend,
                                                                   divergence=divergence,
                                                                   hp_tag=hp_tag_found)
                        if potential_events:
                            if not low_mem:
                                block_edge_events += potential_events
                            else:
                                for res in potential_events:
                                    pickle.dump(res, completed_file)
                    else:
                        j_submitted, w_idx = heapq.heappop(minhq)
                        heapq.heappush(minhq, (j_submitted + len(res.n2n), w_idx))
                        msg_queues[w_idx][1].send(res)
            else:
                # most partitions processed here, dict returned, or None
                res = graph.proc_component(node_to_name, component, read_buffer, infile, G, lower_bound_support,
                                           procs, paired_end, sites_index)
                if res:
                    event_id += 1
                    # Res is a dict {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}
                    if procs == 1:
                        potential_events, event_id = component_job(infile, res, regions, event_id, clip_length,
                                                                   insert_median,
                                                                   insert_stdev,
                                                                   insert_ppf,
                                                                   min_support,
                                                                   lower_bound_support,
                                                                   merge_dist,
                                                                   regions_only,
                                                                   assemble_contigs,
                                                                   rel_diffs=rel_diffs, diffs=diffs, min_size=min_size,
                                                                   max_single_size=max_single_size,
                                                                   sites_index=sites_index,
                                                                   paired_end=paired_end,
                                                                   length_extend=length_extend, divergence=divergence,
                                                                   hp_tag=hp_tag_found)
                        if potential_events:
                            if not low_mem:
                                block_edge_events += potential_events
                            else:
                                for res in potential_events:
                                    pickle.dump(res, completed_file)
                    else:
                        j_submitted, w_idx = heapq.heappop(minhq)
                        heapq.heappush(minhq, (j_submitted + len(res.n2n), w_idx))
                        msg_queues[w_idx][1].send(res)

    if completed_file is not None:
        completed_file.close()
    del G
    del read_buffer
    cmp.clear()
    if low_mem:
        os.remove(component_path)
    gc.collect()
    # #
    if procs > 1 or low_mem:
        for w in msg_queues:
            w[1].send(0)
        for n in consumers:
            n.join()
        event_id = 0
        last_seen_grp_id = None
        new_grp = None
        for p in range(procs):
            jf = open(f"{tdir}/job_{p}.done.pkl", "rb")
            while 1:
                try:
                    res = pickle.load(jf)
                    event_id += 1
                    current_id = res.event_id
                    current_grp = res.grp_id
                    if last_seen_grp_id != current_grp:
                        last_seen_grp_id = current_grp
                        new_grp = event_id  # replace old grp id with new_grp (id of first in group)
                    res.event_id = event_id
                    res.grp_id = new_grp
                    block_edge_events.append(res)
                except EOFError:
                    break
            last_grp_id = event_id
            os.remove(f"{tdir}/job_{p}.done.pkl")
    if len(block_edge_events) == 0:
        return [], None

    logging.info("Number of components: {}. N candidates: {}".format(components_seen, len(block_edge_events)))
    keeps = len([i for i in block_edge_events if i.site_info])
    if keeps:
        logging.info("Number of matching SVs from sites: {}".format(keeps))
    preliminaries = []
    if args["remap"] == "True" and args["contigs"] == "True":
        block_edge_events = re_map.remap_soft_clips(block_edge_events, ref_genome,
                                                    keep_unmapped=True if args["pl"] == "pe" else False,
                                                    min_support=min_support)
        logging.info("Re-alignment of soft-clips done. N candidates: {}".format(len(block_edge_events)))

    block_edge_events = consensus.contig_info(block_edge_events)  # GC info, repetitiveness

    cdef EventResult_t d
    for d in block_edge_events:
        d.type = args["pl"]
        if d.posA > d.posB and d.chrA == d.chrB:
            t = d.posA
            d.posA = d.posB
            d.posB = t

    # Merge across calls
    if args["merge_within"] == "True":
        merged = merge_svs.merge_events(block_edge_events, args["merge_dist"], regions, bool(paired_end), try_rev=False, pick_best=False,
                                        debug=True, min_size=args["min_size"],
                                        max_comparisons=args["max_comparisons"] if "max_comparisons" in args else 100,
                                        procs=args['procs'])
    else:
        merged = block_edge_events
    logging.info("Number of candidate SVs merged: {}".format(len(block_edge_events) - len(merged)))
    logging.info("Number of candidate SVs after merge: {}".format(len(merged)))

    before = len(merged)

    if auto_support:
        if not args["keep_small"]:
            merged = [event for event in merged if (event.svlen >= args["min_size"] or event.chrA != event.chrB) or event.site_info]
    else:
        if not args["keep_small"]:
            merged = [event for event in merged if (event.svlen >= args["min_size"] or event.chrA != event.chrB) and (event.su >= args["min_support"] or event.site_info)]
        else:
            merged = [event for event in merged if (event.su >= args["min_support"] or event.site_info)]

    coverage_analyser = post_call.CoverageAnalyser(tdir)
    preliminaries = coverage_analyser.process_events(merged)

    preliminaries = coverage.get_raw_coverage_information(merged, regions, coverage_analyser, infile, args["max_cov"])

    if auto_support:
        preliminaries = post_call.filter_auto_min_support(preliminaries)
        logging.info("Number of candidate SVs dropped with support < min support: {}".format(before - len(merged)))

    preliminaries = post_call.get_badclip_metric(preliminaries, bad_clip_counter, infile, regions)

    preliminaries = re_map.drop_svs_near_reference_gaps(preliminaries, paired_end, ref_genome, args["drop_gaps"] == "True")
    preliminaries = post_call.ref_repetitiveness(preliminaries, ref_genome)
    preliminaries = post_call.strand_binom_t(preliminaries)

    preliminaries = extra_metrics.find_repeat_expansions(preliminaries, insert_stdev)
    preliminaries = post_call.compressability(preliminaries)
    preliminaries = post_call.get_gt_metric2(preliminaries, args["mode"], True)

    preliminaries = post_call.get_ref_base(preliminaries, ref_genome, args["symbolic_sv_size"])

    preliminaries = extra_metrics.sample_level_density(preliminaries, regions)
    preliminaries = coverage_analyser.normalize_coverage_values(preliminaries)

    n_in_grp = Counter([d.grp_id for d in preliminaries])
    for d in preliminaries:
        d.n_in_grp = n_in_grp[d.grp_id]
    if len(preliminaries) == 0:
        return [], None
    bad_clip_counter.tidy()
    return preliminaries, sites_adder


def cluster_reads(args):
    t0 = time.time()
    np.random.seed(1)
    random.seed(1)
    kind = args["sv_aligns"].split(".")[-1]
    kind = "stdin" if kind == "-" else kind
    opts = {"bam": "rb", "cram": "rc", "sam": "r", "-": "rb", "stdin": "rb"}
    if kind not in opts:
        raise ValueError("Input must be a .bam/cam/sam or stdin")
    if (kind == "stdin" or kind == "sam") and args["regions"] is not None:
        raise ValueError("Cannot use stdin and include (an indexed bam is needed)")

    bam_mode = opts[kind]
    args["bam_mode"] = bam_mode
    pysam.set_verbosity(0)
    infile = pysam.AlignmentFile(args["sv_aligns"], bam_mode, threads=args["procs"],
                                 reference_filename=None if kind != "cram" else args["reference"])
    pysam.set_verbosity(3)
    has_index = True
    try:
        infile.check_index()
    except (ValueError, AttributeError):  # attribute error with sam file
        has_index = False
    if not has_index and args["regions"] is not None:
        logging.info("Input file has no index, but --regions was provided, attempting to index")
        infile.close()
        pysam.index(args["sv_aligns"])
        infile = pysam.AlignmentFile(args["sv_aligns"], bam_mode, threads=args["procs"],
                                     reference_filename=None if kind != "cram" else args["reference"])
    ref_genome = pysam.FastaFile(args["reference"])
    ibam = None
    if args["ibam"] is not None:
        kind2 = args["ibam"].split(".")[-1]
        if kind2 == "stdin" or kind2 == "-" or kind2 not in opts:
            raise ValueError("--ibam must be a .bam/cam/sam file")
        ibam = pysam.AlignmentFile(args["ibam"], opts[kind2], threads=1,
                                   reference_filename=None if kind2 != "cram" else args["reference"])        
                                          
    if "RG" in infile.header:
        rg = infile.header["RG"]
        if "SM" in rg[0]:
            if len(rg) > 1:
                sms = set([i["SM"] for i in rg])
                if len(sms) > 1:
                    logging.warning("Warning: more than one sample in @RG, using first sample (SM) for output: {}".format(rg[0]["SM"]))
            sample_name = rg[0]["SM"]
        else:
            sample_name = os.path.splitext(os.path.basename(args["sv_aligns"]))[0]
            logging.warning("Warning: no SM tag in @RG (read-group), using input file name as sample name for output: {}".format(sample_name))
    else:
        sample_name = os.path.splitext(os.path.basename(args["sv_aligns"]))[0]
        logging.warning("Warning: no @RG, using input file name as sample name for output: {}".format(sample_name))

    logging.info("Sample name: {}".format(sample_name))
    try:
        args["thresholds"] = dict(zip(["DEL", "INS", "INV", "DUP", "TRA"], (map(float, args["thresholds"].split(",")))))
    except:
        raise ValueError("--thresholds parameter not understood, require a comma-separated string for DEL,INS,INV,DUP,TRA")
    if "insert_median" not in args and "I" in args:
        im, istd = map(float, args["I"].split(","))
        args["insert_median"] = im
        args["insert_stdev"] = istd
    if args["svs_out"] == "-" or args["svs_out"] is None:
        logging.info("Writing vcf to stdout")
        outfile = stdout
    else:
        logging.info("Writing SVs to {}".format(args["svs_out"]))
        outfile = open(args["svs_out"], "w")
    logging.info("Running pipeline")
    _debug_k = []
    regions = io_funcs.overlap_regions(args["regions"])

    contig_header_lines = ""
    for item in infile.header["SQ"]:
        contig_header_lines += f"\n##contig=<ID={item['SN']},length={item['LN']}>"

    #####################
    #  Run dysgu here   #
    #####################
    events, site_adder = pipe1(args, infile, kind, regions, ibam, ref_genome, sample_name)
    if not events:
        logging.critical("No events found")
        outfile.write(io_funcs.get_header(contig_header_lines) + "\n")
        return

    df = pd.DataFrame.from_records([to_dict(e) for e in events])
    df = post_call.apply_model(df, args["pl"], args["contigs"], args["diploid"], args["thresholds"])
    if args["sites"]:
        df = post_call.update_prob_at_sites(df, events, args["thresholds"], parse_probs=args["parse_probs"] == "True",
                                        default_prob=args["sites_prob"])
        df["site_id"] = ["." if not s else s.id for s in df["site_info"]]
        if args["all_sites"] == "True":
            df = sites_utils.append_uncalled(df, site_adder, infile, parse_probs=args["parse_probs"] == "True")

    if len(df) > 0:

        df = df.sort_values(["chrA", "posA", "event_id"])
        df["sample"] = [sample_name] * len(df)
        df.rename(columns={"contig": "contigA", "contig2": "contigB"}, inplace=True)
        if args["out_format"] == "csv":
            small_output_f = True
            if args["metrics"]:
                small_output_f = False
            cols = io_funcs.col_names(small_output_f)
            fmt = cols.pop()
            cols += fmt
            df[cols].to_csv(outfile, index=False)
        else:
            args["add_kind"] = "True"
            args["sample_name"] = sample_name
            io_funcs.to_vcf(df, args, {sample_name}, outfile, show_names=False, contig_names=contig_header_lines,
                            sort_output=False)
    else:
        outfile.write(io_funcs.get_header(contig_header_lines) + "\n")

    logging.info("dysgu call {} complete, n={}, time={} h:m:s".format(
               args["sv_aligns"],
               len(df),
               str(datetime.timedelta(seconds=int(time.time() - t0)))))
