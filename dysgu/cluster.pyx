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
from dysgu cimport map_set_utils


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
    # max_dist = max_dist * 1.5
    interval_table = defaultdict(lambda: graph.Table())

    # Make a nested containment list for lookups
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
            if idx != jdx:
                ej = potential[jdx]

                yield ei, ej, idx, jdx


def enumerate_events(G, potential, max_dist, try_rev, tree, rel_diffs=False, diffs=15):


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

        any_ci = bool(ei["contig"] or ei["contig2"])
        any_cj = bool(ej["contig"] or ej["contig2"])

        if loci_similar:

            if any_ci and any_cj and loci_same:  # Try and match contigs

                # Each breakpoint can have a different assembly, only check for match if contigs overlap
                if assembler.check_contig_match(ci, cj, rel_diffs, diffs=diffs, return_int=True) == 1:
                    G.add_edge(idx, jdx)
                    continue

                if assembler.check_contig_match(ci2, cj2, rel_diffs, diffs=diffs, return_int=True) == 1:
                    G.add_edge(idx, jdx)
                    continue

                if assembler.check_contig_match(ci, cj2, rel_diffs, diffs=diffs, return_int=True) == 1:
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


def merge_events(potential, max_dist, tree, try_rev=False, pick_best=False, add_partners=False,
                 rel_diffs=False, diffs=15):
    """Try and merge similar events, use overlap of both breaks points
    """

    if len(potential) <= 1:
        return potential

    G = nx.Graph()
    enumerate_events(G, potential, max_dist, try_rev, tree, rel_diffs, diffs)

    found = []
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item["event_id"]):
            found.append(item)

    for grp in nx.algorithms.components.connected_components(G):

        c = [potential[n] for n in grp]

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

    potential_events = filter_potential(potential_events, regions, regions_only)

    potential_events = merge_events(potential_events, merge_dist, regions,
                                    try_rev=False, pick_best=True, rel_diffs=rel_diffs, diffs=diffs)
    # echo("final potential events", potential_events)
    return potential_events, event_id


def pipe1(args, infile, kind, regions):
    t6 = time.time()
    # inputbam, int max_cov, tree, read_threads, buffer_size
    regions_only = False if args["regions_only"] == "False" else True

    paired_end = int(args["paired"] == "True")
    assemble_contigs = int(args["contigs"] == "True")

    genome_scanner = coverage.GenomeScanner(infile, args["max_cov"], args["include"], args["procs"],
                                            args["buffer_size"], regions_only,
                                            kind == "stdin")
    insert_median, insert_stdev, read_len = -1, -1, -1
    if args["template_size"] != "":
        try:
            insert_median, insert_stdev, read_len = list(map(int, args["template_size"].split(",")))
        except:
            raise ValueError("Template-size must be in the format 'INT,INT,INT'")

    if paired_end:
        insert_median, insert_stdev = genome_scanner.get_read_length(args["max_tlen"], insert_median, insert_stdev,
                                                                     read_len)

        max_dist = insert_median + (insert_stdev * 5)
        max_clust_dist = 1 * (int(insert_median + (5 * insert_stdev)))

        if args["merge_dist"] is None:
            args["merge_dist"] = max_clust_dist

        click.echo(f"Max clustering dist {max_clust_dist}", err=True)

    else:
        max_dist, max_clust_dist = 10000, 10000
        if args["merge_dist"] is None:
            args["merge_dist"] = 500


    event_id = 0
    block_edge_events = []

    echo("begin time", time.time() - t6)

    min_support = args["min_support"]
    read_buffer = genome_scanner.read_buffer



    t5 = time.time()
    G, node_to_name = graph.construct_graph(genome_scanner,
                                            infile,
                                            max_dist=max_clust_dist,
                                            clustering_dist=max_clust_dist,
                                            minimizer_dist=8,
                                            minimizer_support_thresh=args["z_depth"],
                                            minimizer_breadth=args["z_breadth"],
                                            k=21,
                                            m=10,
                                            clip_l=args["clip_length"],
                                            min_support=min_support,
                                            procs=args["procs"],
                                            mapq_thresh=args["mq"],
                                            debug=None,
                                            paired_end=paired_end
                                            )
    echo("graph time", time.time() - t5)

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

    for start_i, end_i in cc:
        cnt += 1
        if end_i - start_i >= min_support:
            component = list(cmp[start_i: end_i])


            res = graph.proc_component(node_to_name, component, read_buffer, infile, G, min_support)
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
    # echo("-----------")
    # if args["merge_within"] == "True":
    #     merged = merge_events(block_edge_events, args["merge_dist"], regions,
    #                             try_rev=False, pick_best=True)
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

    # if args["procs"] > 1:
    #     raise ValueError("Sorry, only single process is supported currently")
    click.echo("Input file: {} ({}). Processes={}".format(args["sv_aligns"], kind, args["procs"]), err=True)

    infile = pysam.AlignmentFile(args["sv_aligns"], opts[kind], threads=args["procs"])

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

    #
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
    events, extended_tags = pipe1(args, infile, kind, regions)
    echo("evet time", time.time() - t4)
    df = pd.DataFrame.from_records(events)



    if len(df) > 0:
        df = df.sort_values(["kind", "chrA", "posA"])

        df["sample"] = [sample_name] * len(df)
        df["id"] = range(len(df))
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
    click.echo("dysgu call {} complete, n={}, mem={} Mb, time={} h:m:s".format(
        args["sv_aligns"],
        len(df),
        int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
        str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)
