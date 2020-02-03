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
from dysgu import coverage, graph, call_component, assembler, io_funcs

def echo(*args):
    click.echo(args, err=True)


def filter_potential(input_events, tree, regions_only):
    potential = []

    for i in input_events:
        if "posB" not in i:  # Skip events for which no posB was identified
            continue
        if i["chrA"] == i["chrB"] and i["posA"] == i["posB"]:
            continue
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
    max_dist = max_dist * 1.5
    interval_table = defaultdict(lambda: graph.Table())

    # Make a nested containment list for faster lookups
    for idx in range(len(potential)):
        ei = potential[idx]
        interval_table[ei["chrA"]].add(ei["posA"] - max_dist, ei["posA"] + max_dist, idx)
        interval_table[ei["chrB"]].add(ei["posB"] - max_dist, ei["posB"] + max_dist, idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}

    for idx in range(len(potential)):
        ei = potential[idx]

        # Overlap right hand side # list(tree[chrom].find_overlap(start, end))
        ols = list(nc[ei["chrB"]].find_overlap(ei["posB"], ei["posB"] + 1))
        for target in ols:
            jdx = target[2]
            ej = potential[jdx]
            yield ei, ej, idx, jdx


def compare_all(potential):
    # Quadratic run time but fast for a hand full of events
    for idx in range(len(potential)):
        ei = potential[idx]

        for jdx in range(len(potential)):
            ej = potential[jdx]

            yield ei, ej, idx, jdx


def enumerate_events(potential, max_dist, try_rev, tree):

    G = nx.Graph()

    if len(potential) < 3:
        event_iter = compare_all(potential)
    else:
        event_iter = compare_subset(potential, max_dist)

    seen = set([])
    for ei, ej, idx, jdx in event_iter:
        i_id = ei["event_id"]
        j_id = ej["event_id"]

        if i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen:
            continue

        seen.add((i_id, j_id))

        # Check if events point to the same loci
        loci_similar = False
        loci_same = False
        if ei["chrA"] == ej["chrA"]:  # Try chrA matches chrA

            dist1 = abs(ei["posA"] - ej["posA"])
            if dist1 < max_dist:
                if ei["chrB"] == ej["chrB"]:

                    dist2 = abs(ei["posB"] - ej["posB"])
                    if dist2 < max_dist:
                        loci_similar = True
                    if dist1 < 5 and dist2 < 5:
                        loci_same = True

        if not loci_similar:  # Try chrA matches chrB
            if ei["chrA"] == ej["chrB"]:
                dist1 = abs(ei["posA"] - ej["posB"])
                if dist1 < max_dist:
                    if ei["chrB"] == ej["chrA"]:
                        dist2 = abs(ei["posB"] - ej["posA"])
                        if dist2 < max_dist:
                            loci_similar = True
                        if dist1 < 5 and dist2 < 5:
                            loci_same = True

        if "contig" not in ei or "contig" not in ej:
            continue

        ci = "" if ei["contig"] is None else ei["contig"]
        ci2 = "" if ei["contig2"] is None else ei["contig2"]
        cj = "" if ej["contig"] is None else ej["contig"]
        cj2 = "" if ej["contig2"] is None else ej["contig2"]

        any_ci = ei["contig"] or ei["contig2"]
        any_cj = ej["contig"] or ej["contig2"]


        if loci_similar:

            if any_ci and any_cj and loci_same:  # Try and match contigs

                # Each breakpoint can have a different assembly, only check for match if contigs overlap
                if assembler.check_contig_match(ci, cj, diffs=15, return_int=True) == 1:
                    G.add_edge(idx, jdx)
                    continue

                if assembler.check_contig_match(ci2, cj2, diffs=15, return_int=True) == 1:
                    G.add_edge(idx, jdx)
                    continue

                if assembler.check_contig_match(ci, cj2, diffs=15, return_int=True) == 1:
                    G.add_edge(idx, jdx)
                    continue

            # No contigs to match, merge anyway
            elif not any_ci and not any_cj and loci_same:
                G.add_edge(idx, jdx)  #G.add_edge(i_id, j_id)

            # Only merge loci if they are not both within --include regions. chrA:posA only needs checking
            elif not (io_funcs.intersecter_str_chrom(tree, ei["chrA"], ei["posA"], ei["posA"] + 1) and
                      io_funcs.intersecter_str_chrom(tree, ei["chrB"], ei["posB"], ei["posB"] + 1)):
                G.add_edge(idx, jdx)
    return G


def merge_events(potential, max_dist, tree, try_rev=False, pick_best=False, add_partners=False):

    if len(potential) <= 1:
        return potential  #, seen

    G = enumerate_events(potential, max_dist, try_rev, tree)

    found = []
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item["event_id"]):
            found.append(item)

    for grp in nx.algorithms.components.connected_components(G):
        try:
            c = [potential[n] for n in grp]
        except:
            echo("WARN", len(grp), len(potential))
            continue
            # echo(grp.nodes())
            # quit()

        best = sorted(c, key=lambda x: sum([x["pe"], x["supp"]]), reverse=True)
        w0 = best[0]["pe"] + best[0]["supp"]  # Weighting for base result

        if not pick_best:
            for k in range(1, len(best)):

                # Sum these
                for t in ["pe", "supp", "sc", "NP", "block_edge", "joinA", "joinB"]:
                    best[0][t] += best[k][t]

                if best[k]["maxASsupp"] > best[0]["maxASsupp"]:
                    best[0]["maxASsupp"] = best[k]["maxASsupp"]

                # Average these
                for t in ["DN", "MAPQsupp", "MAPQpri", "DApri", "DAsupp", "DP", "NMpri", "NMsupp"]:
                    w = best[k]["pe"] + best[k]["supp"]
                    denom = w0 + w
                    if denom == 0:
                        weighted_av = 0
                    else:
                        weighted_av = ((best[0][t] * w0) + (best[k][t] * w)) / denom
                    best[0][t] = weighted_av
                w0 = best[0]["pe"] + best[0]["supp"]

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
        if not io_funcs.intersecter_int_chrom(regions, ei["posA"], ei["posA"], ei["posA"] + 1):
            interval_table[ei["chrA"]].add(ei["posA"] - max_dist, ei["posA"] + max_dist, idx)

        if not io_funcs.intersecter_int_chrom(regions, ei["posB"], ei["posB"], ei["posB"] + 1):
            interval_table[ei["chrB"]].add(ei["posB"] - max_dist, ei["posB"] + max_dist, idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}

    for idx in range(len(potential)):
        ei = potential[idx]

        # Overlap right hand side # list(tree[chrom].find_overlap(start, end))
        neighbors = 0.
        count = 0.
        if not io_funcs.intersecter_int_chrom(regions, ei["posA"], ei["posA"], ei["posA"] + 1):
            neighbors += len(list(nc[ei["chrA"]].find_overlap(ei["posA"], ei["posA"] + 1))) - 1
            count += 1

        if not io_funcs.intersecter_int_chrom(regions, ei["posB"], ei["posB"], ei["posB"] + 1):
            neighbors += len(list(nc[ei["chrB"]].find_overlap(ei["posB"], ei["posB"] + 1))) - 1
            count += 1

        if neighbors < 0:
            neighbors = 0
        if count > 0:
            ei["neigh"] = neighbors / count
        else:
            ei["neigh"] = 0
    return potential


def component_job(infile, component, regions, event_id, max_dist, clip_length, insert_med, insert_stdev, min_supp,
                  regions_only):

    potential_events = []

    for event in call_component.call_from_block_model(infile,
                                                      component,
                                                      clip_length,
                                                      insert_med,
                                                      insert_stdev,
                                                      min_supp):
        if event:
            event["event_id"] = event_id
            if event["chrA"] is not None:
                event["chrA"] = infile.get_reference_name(event["chrA"])
            if event["chrB"] is not None:
                event["chrB"] = infile.get_reference_name(event["chrB"])

            potential_events.append(event)

    potential_events = filter_potential(potential_events, regions, regions_only)

    potential_events = merge_events(potential_events, max_dist, regions,
                                    try_rev=False, pick_best=True)

    return potential_events


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory {}".format(d))


def cluster_reads(args):
    t0 = time.time()
    np.random.seed(1)
    random.seed(1)
    try:
        model = pickle.load(open(args["model"], "rb"))
        click.echo("Model loaded from {}".format(args["model"]), err=True)
    except UserWarning:
        model = pickle.load(open(args["model"], "rb"))
        click.echo("Model loaded, outdated model vs sklearn? {}".format(args["model"]), err=True)
    except:
        model = None
        click.echo("No model loaded", err=True)

    mk_dest(args["dest"])
    if args["dest"] is None:
        args["dest"] = "."

    kind = args["sv_aligns"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "rs"}

    click.echo("Input file is {}, (.{} format). Processes={}".format(args["sv_aligns"], kind, args["procs"]), err=True)

    infile = pysam.AlignmentFile(args["sv_aligns"], opts[kind], threads=args["procs"])

    sample_name = os.path.splitext(os.path.basename(args["sv_aligns"]))[0]

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

    def pipe1():

        t0 = time.time()
        # inputbam, int max_cov, tree, read_threads, buffer_size
        regions_only = False if args["regions_only"] == "False" else True

        genome_scanner = coverage.GenomeScanner(infile, args["max_cov"], args["include"], args["procs"],
                                                   args["buffer_size"], regions_only)
        insert_median, insert_stdev, read_len = -1, -1, -1
        if args["template_size"] != "":
            try:
                insert_median, insert_stdev, read_len = list(map(int, args["template_size"].split(",")))
            except:
                raise ValueError("Template-size must be in the format 'INT,INT,INT'")

        insert_median, insert_stdev = genome_scanner.get_read_length(args["max_tlen"], insert_median, insert_stdev,
                                                                     read_len)


        max_clust_dist = 1 * (int(insert_median + (5 * insert_stdev)))  # 2
        click.echo(f"Max clustering dist {max_clust_dist}", err=True)
        event_id = 0
        block_edge_events = []

        for component in graph.construct_graph(genome_scanner,
                                                  infile,
                                                  max_dist=insert_median + (insert_stdev * 5),
                                                  clustering_dist=max_clust_dist,
                                                  minimizer_dist=8,
                                                  minimizer_support_thresh=args["z_depth"],
                                                  minimizer_breadth=args["z_breadth"],
                                                  k=21,
                                                  m=10,
                                                  clip_l=args["clip_length"],
                                                  min_support=args["min_support"],
                                                  procs=args["procs"],
                                                  debug=None):

            potential_events = component_job(infile, component, regions, event_id, max_clust_dist, args["clip_length"],
                                                      insert_median,
                                                      insert_stdev,
                                                      args["min_support"],
                                                      regions_only)

            event_id += 1
            block_edge_events += potential_events

        # Need to parallelize after the reads are collected, otherwise multiple file handles causes problems
        # else:
            # fileinfo = (args["sv_aligns"], opts[kind])
            #
            # with parallel_backend("multiprocessing"):
            #
            #     block_edges = Parallel(n_jobs=args["procs"],
            #                            )(delayed(component_job)(fileinfo,  #infile,
            #                                                      component,
            #                                                      regions,
            #                                                      event_id,
            #                                                      max_dist,
            #                                                      args["clip_length"],
            #                                                      args["insert_median"],
            #                                                      args["insert_stdev"],
            #                                                      args["min_support"])
            #                               for event_id, component in enumerate(cy_graph.construct_graph(infile,
            #                                               int(args["insert_median"]),
            #                                               int(max_dist),
            #                                               regions,
            #                                               minimizer_dist=args["read_length"],
            #                                               k=21,
            #                                               m=8,
            #                                               clip_l=args["clip_length"],
            #                                               input_windows=regions_merged,
            #                                               buf_size=args["buffer_size"],
            #                                               min_support=args["min_support"],
            #                                               procs=args["procs"],
            #                                               debug=None)))
            #
            # block_edge_events = [item for sublist in block_edges for item in sublist]

        click.echo("Processed chunks {}s".format(round(time.time() - t0, 1)), err=True)

        t0 = time.time()

        preliminaries = []

        # Merge across calls
        # merged = merge_events(block_edge_events, max_clust_dist / 2., regions,
        #                             try_rev=False, pick_best=True)
        merged = block_edge_events
        if merged:
            for event in merged:
                # Collect coverage information
                event_dict = call_component.get_raw_coverage_information(event, regions, genome_scanner.depth_d, infile) # regions_depth)

                if event_dict:
                    preliminaries.append(event_dict)


        preliminaries = assembler.contig_info(preliminaries)  # GC info, pretetitiveness
        preliminaries = sample_level_density(preliminaries, regions)

        return preliminaries

    preliminaries = pipe1()

    classified_events_df = call_component.calculate_prob_from_model(preliminaries, model)

    # "mark", "mark_seq", "mark_ed", "templated_ins_info",
    #          "templated_ins_len",
    # Out order
    k = ["chrA", "posA", "chrB", "posB", "sample", "id", "kind", "svtype", "join_type", "cipos95A", "cipos95B",
         "DP", "DN", "DApri", "DAsupp",  "NMpri", "NMsupp", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "pe", "supp", "sc", "block_edge",
         "raw_reads_10kb",
          "linked", "contigA", "contigB",  "gc", "neigh", "rep", "ref_bases", "Prob"]

    c = 0
    if classified_events_df is not None and len(classified_events_df) > 0:
        c = len(classified_events_df)
        classified_events_df["sample"] = [sample_name]*len(classified_events_df)
        classified_events_df["id"] = range(len(classified_events_df))
        classified_events_df = classified_events_df.rename(columns={"contig": "contigA", "contig2": "contigB"})
        classified_events_df[k].to_csv(outfile, index=False)

    click.echo("call-events {} complete, n={}, {} h:m:s, {}".format(args["sv_aligns"],
                                                                c,
                                                                str(datetime.timedelta(seconds=int(time.time() - t0))),
                                                                time.time() - t0),
               err=True)
