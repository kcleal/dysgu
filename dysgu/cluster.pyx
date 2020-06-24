# cython: language_level=3

from __future__ import absolute_import
import datetime
import os
import time
import numpy as np
import random
from collections import defaultdict
import click
import networkx as nx
import pysam
import sys
import pickle
import resource
import pandas as pd
from dysgu import coverage, graph, call_component, assembler, io_funcs
from dysgu.map_set_utils cimport is_reciprocal_overlapping, span_position_distance, position_distance
import itertools


def echo(*args):
    click.echo(args, err=True)


def filter_potential(input_events, tree, regions_only):
    potential = []

    for i in input_events:

        if "posB" not in i:  # Skip events for which no posB was identified
            continue
        # if i["chrA"] == i["chrB"] and i["posA"] == i["posB"]:
        #     continue
        if "contig" not in i or i["contig"] == "":
            i["contig"] = None

        # Remove events for which both ends are in --include but no contig was found

        posA_intersects = io_funcs.intersecter_str_chrom(tree, i["chrA"], i["posA"], i["posA"] + 1)
        posB_intersects = io_funcs.intersecter_str_chrom(tree, i["chrB"], i["posB"], i["posB"] + 1)
        if (posA_intersects and posB_intersects) and (i["contig"] is None or i["contig2"] is None):
            continue

        # Remove events for which neither end is in --include (if --include provided)
        if tree and regions_only:
            if not posA_intersects and not posB_intersects:
                continue

        potential.append(i)
    return potential


def compare_subset(potential, max_dist):

    interval_table = defaultdict(lambda: graph.Table())

    # Make a nested containment list for lookups
    for idx in range(len(potential)):
        ei = potential[idx]

        interval_table[ei["chrA"]].add(ei["posA"] - ei["cipos95A"] - max_dist, ei["posA"] + ei["cipos95A"] + max_dist, idx)
        if ei["svtype"] != "INS":
            if ei["chrA"] != ei["chrB"]:
                interval_table[ei["chrB"]].add(ei["posB"] - ei["cipos95B"] - max_dist, ei["posB"] + ei["cipos95B"] + max_dist, idx)
            else:
                dist1 = abs(ei["posB"] - ei["posA"])
                ci_a = max(ei["cipos95A"], ei["cipos95A"])
                if dist1 + ci_a > max_dist:
                    interval_table[ei["chrB"]].add(ei["posB"] - ei["cipos95B"] - max_dist, ei["posB"] + ei["cipos95B"] + max_dist, idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}

    for idx in range(len(potential)):
        ei = potential[idx]

        # Overlap right hand side
        ols = list(nc[ei["chrB"]].find_overlap(ei["posB"], ei["posB"] + 1))
        for target in ols:
            jdx = target[2]
            ej = potential[jdx]
            yield ei, ej, idx, jdx


def compare_all(potential):
    # Quadratic run time but fast for a hand full of events
    for idx, jdx in itertools.product(range(len(potential)), range(len(potential))):
        yield potential[idx], potential[jdx], idx, jdx


def enumerate_events(G, potential, max_dist, try_rev, tree, rel_diffs=False, diffs=15, same_sample=True, debug=False):

    if len(potential) < 3:
        event_iter = compare_all(potential)  # N^2 compare all events to all others
        # echo("compare all")
    else:
        event_iter = compare_subset(potential, max_dist)  # Use NCLS, generate overlap tree and perform intersections

    seen = set([])

    for ei, ej, idx, jdx in event_iter:

        i_id = ei["event_id"]
        j_id = ej["event_id"]

        if i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen:
            continue

        seen.add((i_id, j_id))

        if ei["svtype"] != ej["svtype"]:
            continue

        if not same_sample and "sample" in ei:
            if ei["sample"] == ej["sample"]:
                continue

        # Check if events point to the same loci
        loci_similar = False
        loci_same = False

        if ei["chrA"] == ej["chrA"] and ei["chrB"] == ej["chrB"]:  # Try chrA matches chrA

            # if ei["svtype"] != "INS":
                if is_reciprocal_overlapping(ei["posA"] - ei["cipos95A"] - 5,
                                             ei["posB"] + ei["cipos95B"] + 5,
                                             ej["posA"] - ej["cipos95A"] - 5,
                                             ej["posB"] + ej["cipos95B"] + 5):
                    loci_similar = True

                dist1 = abs(ei["posA"] - ej["posA"])
                dist2 = abs(ei["posB"] - ej["posB"])
                if dist1 < 5 and dist2 < 5:
                    loci_same = True

            # else:
            #     dist1 = abs(ei["posA"] - ej["posA"])
            #     dist2 = abs(ei["posB"] - ej["posB"])
            #     echo(ei, ej)
            #     if dist1 < 5 and dist2 < 5:
            #         loci_same = True

            # ci_a = max(ei["cipos95A"], ej["cipos95A"])
            # if dist1 < max_dist:
            #     if ei["chrB"] == ej["chrB"]:
            #         dist2 = abs(ei["posB"] - ej["posB"])
            #         # ci_b = max(ei["cipos95B"], ej["cipos95B"])
            #         # if dist2 < (max_dist + ci_b):
            #         if ei["svtype"] != "INS":
            #             echo(is_reciprocal_overlapping(ei["posA"], ei["posB"], ej["posA"], ej["posB"]))
            #             if is_reciprocal_overlapping(ei["posA"], ei["posB"], ej["posA"], ej["posB"]):
            #                 loci_similar = True
            #         if dist1 < 5 and dist2 < 5:
            #             loci_same = True

        if not loci_similar:  # Try chrA matches chrB

            if ei["chrA"] == ej["chrB"] and ei["chrB"] == ej["chrA"]:  # Try chrA matches chrA

                if is_reciprocal_overlapping(ei["posA"] - ei["cipos95A"] - 5,
                                             ei["posB"] + ei["cipos95B"] + 5,
                                             ej["posA"] - ej["cipos95A"] - 5,
                                             ej["posB"] + ej["cipos95B"] + 5):
                    loci_similar = True

                dist1 = abs(ei["posA"] - ej["posB"])
                dist2 = abs(ei["posB"] - ej["posA"])
                if dist1 < 5 and dist2 < 5:
                    loci_same = True

                # if ei["svtype"] != "INS":
                #     if is_reciprocal_overlapping(ei["posA"], ei["posB"], ej["posA"], ej["posB"]):
                #         loci_similar = True
                # else:
                #     dist1 = abs(ei["posA"] - ej["posB"])
                #     dist2 = abs(ei["posB"] - ej["posA"])
                #     if dist1 < 5 and dist2 < 5:
                #         loci_same = True
                #         loci_similar = True


            # if ei["chrA"] == ej["chrB"]:
            #
            #
            #     dist1 = abs(ei["posA"] - ej["posB"])
            #     # ci_a = max(ei["cipos95A"], ej["cipos95B"])
            #     if dist1 < max_dist:
            #         if ei["chrB"] == ej["chrA"]:
            #             dist2 = abs(ei["posB"] - ej["posA"])
            #             # ci_b = max(ei["cipos95B"], ej["cipos95A"])
            #             if ei["svtype"] != "INS":
            #                 if is_reciprocal_overlapping(ei["posA"], ei["posB"], ej["posA"], ej["posB"]):
            #                     loci_similar = True
            #
            #             # if dist2 < (max_dist + ci_b):
            #             #     loci_similar = True
            #             if dist1 < 5 and dist2 < 5:
            #                 loci_same = True


        if "contig" not in ei or "contig" not in ej:
            continue

        ci = "" if ei["contig"] is None else ei["contig"]
        ci2 = "" if ei["contig2"] is None else ei["contig2"]
        cj = "" if ej["contig"] is None else ej["contig"]
        cj2 = "" if ej["contig2"] is None else ej["contig2"]

        if ei["svtype"] != "INS":
            any_contigs_to_check = ei["contig"] or ei["contig2"] and (ej["contig"] or ej["contig2"])
        else:
            # No guarentee contigs or in reference order
            any_contigs_to_check = ei["contig"] and ej["contig"] and ((ei["contig2"] and ej["contig2"]) or (not ei["contig2"] and not ej["contig2"]))

        # echo(ei["posA"], ej["posA"], loci_same, loci_similar, assembler.check_contig_match(ci, cj, rel_diffs, diffs=diffs, return_int=True),
        #      any_contigs_to_check)
        if loci_similar or loci_same:

            if any_contigs_to_check and loci_same:  # Try and match contigs

                # Each breakpoint can have a different assembly, only check for match if contigs overlap
                if assembler.check_contig_match(ci, cj, rel_diffs, diffs=diffs, return_int=True) == 1:
                    G.add_edge(i_id, j_id, loci_same=loci_same)
                    continue

                if assembler.check_contig_match(ci2, cj2, rel_diffs, diffs=diffs, return_int=True) == 1:
                    G.add_edge(i_id, j_id, loci_same=loci_same)
                    continue

                if assembler.check_contig_match(ci, cj2, rel_diffs, diffs=diffs, return_int=True) == 1:
                    G.add_edge(i_id, j_id, loci_same=loci_same)
                    continue

            # No contigs to match, merge anyway
            elif not any_contigs_to_check and loci_same:
                G.add_edge(i_id, j_id, loci_same=loci_same)

            # Only merge loci if they are not both within --include regions. chrA:posA only needs checking
            elif not (io_funcs.intersecter_str_chrom(tree, ei["chrA"], ei["posA"], ei["posA"] + 1) and
                      io_funcs.intersecter_str_chrom(tree, ei["chrB"], ei["posB"], ei["posB"] + 1)):
                G.add_edge(i_id, j_id, loci_same=loci_same)
    return G


def cut_components(G):
    e = G.edges(data=True)
    G2 = nx.Graph([i for i in e if i[2]["loci_same"] == True])
    for u in G.nodes():
        if u not in G2:
            e0 = next(G.edges(u).__iter__())  # Use first edge out of u to connect
            G2.add_edge(*e0)
    return nx.algorithms.components.connected_components(G2)


cpdef srt_func(c):
    if "su" in c:
        return c["su"]
    # over weight spanning
    return c["pe"] + c["supp"] + (c["spanning"] * 10) + c["sc"]


def merge_events(potential, max_dist, tree, try_rev=False, pick_best=False, add_partners=False,
                 rel_diffs=False, diffs=15, same_sample=True, debug=False):
    """Try and merge similar events, use overlap of both breaks points
    """
    max_dist = max_dist / 2
    if len(potential) <= 1:
        return potential

    G = nx.Graph()
    G = enumerate_events(G, potential, max_dist, try_rev, tree, rel_diffs, diffs, same_sample, debug=debug)  # Cluster events on graph

    found = []
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item["event_id"]):
            found.append(item)

    # discard = set([])  # event_id's that were merged, delete these later

    # Try and merge SVs with identical breaks, then merge ones with less accurate breaks - this helps prevent
    # over merging SVs that are close together

    components = cut_components(G)
    node_to_event = {i["event_id"]: i for i in potential}
    cdef int k
    # Only keep edges with loci_same==False if removing the edge leads to an isolated node
    for grp in components:

        c = [node_to_event[n] for n in grp]
        best = sorted(c, key=srt_func, reverse=True)

        if not pick_best:

            w0 = best[0]  # Weighting for base result
            weight = w0["pe"] + w0["supp"] + w0["spanning"]
            for k in range(1, len(best)):

                item = best[k]

                # Sum these
                for t in ("pe", "supp", "sc", "su", "NP", "block_edge", "plus", "minus", "spanning"):
                    if t in item:
                        best[0][t] += item[t]

                if item["maxASsupp"] > w0["maxASsupp"]:
                    best[0]["maxASsupp"] = item["maxASsupp"]

                # Average weight these
                wt = item["pe"] + item["supp"] + item["spanning"]
                for t in ("DN", "MAPQsupp", "MAPQpri", "DApri", "DAsupp", "DP", "NMpri", "NMsupp"):
                    if t in best[0]:
                        denom = weight + wt
                        if denom == 0:
                            weighted_av = 0
                        else:
                            weighted_av = ((w0[t] * weight) + (item[t] * wt)) / denom
                        best[0][t] = weighted_av

        if add_partners:
            best[0]["partners"] = [i["event_id"] for i in best[1:]]

        found.append(best[0])

    return found


def sample_level_density(potential, regions, max_dist=50):

    interval_table = defaultdict(lambda: graph.Table())

    # Make a nested containment list for faster lookups
    for idx in range(len(potential)):
        ei = potential[idx]

        # Only find density for non-region calls, otherwise too dense to be meaningful
        if not io_funcs.intersecter_str_chrom(regions, ei["chrA"], ei["posA"], ei["posA"] + 1):
            interval_table[ei["chrA"]].add(ei["posA"] - max_dist, ei["posA"] + max_dist, idx)

        if not io_funcs.intersecter_str_chrom(regions, ei["chrB"], ei["posB"], ei["posB"] + 1):
            interval_table[ei["chrB"]].add(ei["posB"] - max_dist, ei["posB"] + max_dist, idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}

    for idx in range(len(potential)):
        ei = potential[idx]

        # Overlap right hand side # list(tree[chrom].find_overlap(start, end))
        neighbors = 0.
        count = 0.
        if ei["chrA"] == ei["chrB"] and abs(ei["posB"] - ei["posA"]) < 2:
            expected = 2
        else:
            expected = 1
        if not io_funcs.intersecter_str_chrom(regions, ei["chrA"], ei["posA"], ei["posA"] + 1):
            neighbors += len(list(nc[ei["chrA"]].find_overlap(ei["posA"], ei["posA"] + 1))) - expected
            count += 1

        if not io_funcs.intersecter_str_chrom(regions, ei["chrB"], ei["posB"], ei["posB"] + 1):
            neighbors += len(list(nc[ei["chrB"]].find_overlap(ei["posB"], ei["posB"] + 1))) - expected
            count += 1

        # Merge overlapping intervals
        neighbors_10kb = 0.
        count_10kb = 0
        large_itv = coverage.merge_intervals(((ei["chrA"], ei["posA"], ei["posA"] + 1),
                                              (ei["chrB"], ei["posB"], ei["posB"] + 1)), pad=10000)

        for c, s, e in large_itv:
            if not io_funcs.intersecter_str_chrom(regions, c, s, e):
                neighbors_10kb += len(list(nc[c].find_overlap(s, e))) - len(large_itv)
                count_10kb += 1

        if neighbors < 0:
            neighbors = 0
        if count > 0:
            ei["neigh"] = neighbors / count
        else:
            ei["neigh"] = 0
        if count_10kb > 0:
            ei["neigh10kb"] = neighbors_10kb / count_10kb
        else:
            ei["neigh10kb"] = 0

    return potential


def component_job(infile, component, regions, event_id, max_dist, clip_length, insert_med, insert_stdev, min_supp,
                  merge_dist, regions_only, extended_tags, assemble_contigs, rel_diffs, diffs):

    potential_events = []
    for event in call_component.call_from_block_model(infile,
                                                      component,
                                                      clip_length,
                                                      insert_med,
                                                      insert_stdev,
                                                      min_supp,
                                                      extended_tags,
                                                      assemble_contigs):
        # echo(event["chrA"], event["posA"])

        if event:
            event["event_id"] = event_id
            if event["chrA"] is not None:
                event["chrA"] = infile.get_reference_name(event["chrA"])
            if event["chrB"] is not None:
                event["chrB"] = infile.get_reference_name(event["chrB"])
            potential_events.append(event)
            event_id += 1
    if not potential_events:
        return potential_events, event_id
    # echo(potential_events)
    potential_events = filter_potential(potential_events, regions, regions_only)

    potential_events = merge_events(potential_events, merge_dist, regions,
                                    try_rev=False, pick_best=True, rel_diffs=rel_diffs, diffs=diffs)

    return potential_events, event_id


def pipe1(args, infile, kind, regions, ibam):

    regions_only = False if args["regions_only"] == "False" else True
    paired_end = int(args["paired"] == "True")
    assemble_contigs = int(args["contigs"] == "True")

    genome_scanner = coverage.GenomeScanner(infile, args["max_cov"], args["include"], args["procs"],
                                            args["buffer_size"], regions_only,
                                            kind == "stdin",
                                            clip_length=args["clip_length"],
                                            min_within_size=args["min_size"])
    insert_median, insert_stdev, read_len = -1, -1, -1
    if args["template_size"] != "":
        try:
            insert_median, insert_stdev, read_len = list(map(int, args["template_size"].split(",")))
        except:
            raise ValueError("Template-size must be in the format 'INT,INT,INT'")

    if paired_end:
        insert_median, insert_stdev = genome_scanner.get_read_length(args["max_tlen"], insert_median, insert_stdev,
                                                                     read_len, ibam)

        max_dist = int(insert_median + (insert_stdev * 5))  # 5
        max_clust_dist = 5 * (int(insert_median + (5 * insert_stdev)))
        if args["merge_dist"] is None:
            args["merge_dist"] = max_clust_dist
        click.echo(f"Max clustering dist {max_clust_dist}", err=True)

    else:
        max_dist, max_clust_dist = 35, 500000
        if args["merge_dist"] is None:
            args["merge_dist"] = 50

    event_id = 0
    block_edge_events = []
    min_support = args["min_support"]
    read_buffer = genome_scanner.read_buffer

    t5 = time.time()
    G, node_to_name = graph.construct_graph(genome_scanner,
                                            infile,
                                            max_dist=max_dist,
                                            clustering_dist=max_clust_dist,
                                            minimizer_dist=8,
                                            minimizer_support_thresh=args["z_depth"],
                                            minimizer_breadth=args["z_breadth"],
                                            k=21,
                                            m=10,
                                            clip_l=args["clip_length"],
                                            min_sv_size=args["min_size"],
                                            min_support=min_support,
                                            procs=args["procs"],
                                            mapq_thresh=args["mq"],
                                            debug=None,
                                            paired_end=paired_end
                                            )
    echo("graph time", time.time() - t5)
    click.echo("graph time, mem={} Mb, time={} h:m:s".format(
        int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
        time.time() - t5), err=True)

    cdef int cnt = 0

    t0 = time.time()
    cmp, cc = graph.get_block_components(G, node_to_name, infile, read_buffer, min_support)

    extended_tags = genome_scanner.extended_tags

    if paired_end:
        rel_diffs = False
        diffs = 15
    else:
        rel_diffs = True
        diffs = 0.15
    echo("len components", len(cc))
    for start_i, end_i in cc:
        cnt += 1
        if end_i - start_i >= min_support:
            component = list(cmp[start_i: end_i])

            res = graph.proc_component(node_to_name, component, read_buffer, infile, G, min_support)
            if res:
                # Res is a dict
                # {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}
                potential_events, event_id = component_job(infile, res, regions, event_id, max_clust_dist, args["clip_length"],
                                                 insert_median,
                                                 insert_stdev,
                                                 args["min_support"],
                                                 args["merge_dist"],
                                                 regions_only,
                                                 extended_tags,
                                                 assemble_contigs,
                                                 rel_diffs=rel_diffs, diffs=diffs)
                if potential_events:
                    block_edge_events += potential_events

    echo("components", time.time() - t0)
    t3 = time.time()
    preliminaries = []

    # Merge across calls
    if args["merge_within"] == "True":
        tmr = time.time()
        block_edge_events = merge_events(block_edge_events, args["merge_dist"], regions, try_rev=False, pick_best=False,
                                         debug=True)
        echo("merging all", time.time() - tmr)
    merged = block_edge_events
    if merged:
        for event in merged:
            # Collect coverage information
            event_dict = call_component.get_raw_coverage_information(event, regions, genome_scanner.depth_d, infile) # regions_depth)

            if event_dict:
                preliminaries.append(event_dict)

    preliminaries = assembler.contig_info(preliminaries)  # GC info, repetitiveness
    preliminaries = sample_level_density(preliminaries, regions)
    echo("time3", time.time() - t3)
    return preliminaries, extended_tags


def cluster_reads(args):
    t0 = time.time()
    np.random.seed(1)
    random.seed(1)
    kind = args["sv_aligns"].split(".")[-1]
    kind = "stdin" if kind == "-" else kind
    opts = {"bam": "rb", "cram": "rc", "sam": "rs", "-": "rb", "stdin": "rb"}
    if kind not in opts:
        raise ValueError("Input must be a .bam/cam/sam or stdin")
    # if args["procs"] > 1:
    #     raise ValueError("Sorry, only single process is supported currently")
    click.echo("Input file: {} ({}). Processes={}".format(args["sv_aligns"], kind, args["procs"]), err=True)
    infile = pysam.AlignmentFile(args["sv_aligns"], opts[kind], threads=args["procs"])

    ibam = None
    if args["ibam"] is not None:
        kind = args["sv_aligns"].split(".")[-1]
        if kind == "stdin" or kind == "-" or kind not in opts:
            raise ValueError("--ibam must be a .bam/cam/sam file")
        ibam = pysam.AlignmentFile(args["ibam"], opts[kind], threads=args["procs"])

    if "RG" in infile.header:
        rg = infile.header["RG"]
        if len(rg) > 1:
            click.echo("Warning: more than one @RG, using first sample (SM) for output: {}".format(rg[0]["SM"]),
                       err=True)
        sample_name = rg[0]["SM"]

    else:
        sample_name = os.path.splitext(os.path.basename(args["sv_aligns"]))[0]
        click.echo("Warning: no @RG, using input file name as sample name for output: {}".format(sample_name), err=True)

    click.echo("Sample name: {}".format(sample_name), err=True)

    if "insert_median" not in args and "I" in args:
        im, istd = map(float, args["I"].split(","))
        args["insert_median"] = im
        args["insert_stdev"] = istd

    if args["svs_out"] == "-" or args["svs_out"] is None:
        click.echo("SVs output to stdout", err=True)
        outfile = sys.stdout
    else:
        click.echo("SVs output to {}".format(args["svs_out"]), err=True)
        outfile = open(args["svs_out"], "w")

    _debug_k = []
    regions = io_funcs.overlap_regions(args["include"])

    # Run dysgu here:
    t4 = time.time()
    events, extended_tags = pipe1(args, infile, kind, regions, ibam)
    echo("evet time", time.time() - t4)
    df = pd.DataFrame.from_records(events)

    if len(df) > 0:
        df = df.sort_values(["kind", "chrA", "posA"])
        df["sample"] = [sample_name] * len(df)
        df["id"] = range(len(df))
        m = {"pe": 1, "pacbio": 2, "nanopore": 3}
        df["type"] = [m[args["mode"]]] * len(df)
        df.rename(columns={"contig": "contigA", "contig2": "contigB"}, inplace=True)
        if args["out_format"] == "csv":
            df[io_funcs.col_names()].to_csv(outfile, index=False)
        else:
            contig_header_lines = ""
            for item in infile.header["SQ"]:
                contig_header_lines += f"\n##contig=<ID={item['SN']},length={item['LN']}>"
            args["add_kind"] = "True"
            args["sample_name"] = sample_name
            io_funcs.to_vcf(df, args, {sample_name}, outfile, show_names=False, contig_names=contig_header_lines,
                            extended_tags=extended_tags)

    infile.close()
    echo(time.time() - t0)
    click.echo("dysgu call {} complete, n={}, time={} h:m:s".format(
               args["sv_aligns"],
               len(df),
               str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)
