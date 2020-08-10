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
from dysgu import coverage, graph, call_component, assembler, io_funcs, re_map
from dysgu.map_set_utils cimport is_reciprocal_overlapping, span_position_distance, position_distance, is_overlapping
import itertools
from operator import mul


def echo(*args):
    click.echo(args, err=True)


def filter_potential(input_events, tree, regions_only):
    potential = []

    for i in input_events:

        if "posB" not in i:  # Skip events for which no posB was identified
            continue

        if i["svtype"] == "INS" and "svlen_precise" in i and i["svlen_precise"] == 0 and not i["contig"] and not i["contig2"]:
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

    # Make a nested containment list for interval intersection
    for idx in range(len(potential)):
        ei = potential[idx]
        interval_table[ei["chrA"]].add(ei["posA"] - ei["cipos95A"] - max_dist, ei["posA"] + ei["cipos95A"] + max_dist, idx)
        # Add another interval for large events, or translocations
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


def enumerate_events(G, potential, max_dist, try_rev, tree, rel_diffs=False, diffs=15, same_sample=True,
                     debug=False):

    if len(potential) < 3:
        event_iter = compare_all(potential)  # N^2 compare all events to all others
        # echo("compare all")
    else:
        event_iter = compare_subset(potential, max_dist)  # Use NCLS, generate overlap tree and perform intersections

    seen = set([])
    pad = 100
    for ei, ej, idx, jdx in event_iter:

        i_id = ei["event_id"]
        j_id = ej["event_id"]
        if i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen:
            continue
        # echo(idx, jdx)
        seen.add((i_id, j_id))

        fail = False
        if ei["svtype"] != ej["svtype"]:
            fail = True
        if fail:
            continue

        if not same_sample and "sample" in ei:
            if ei["sample"] == ej["sample"]:
                continue

        # Check if events point to the same loci
        loci_similar = False
        loci_same = False

        if ei["chrA"] == ej["chrA"] and ei["chrB"] == ej["chrB"]:  # Try chrA matches chrA
            dist1 = abs(ei["posA"] - ej["posA"])
            dist2 = abs(ei["posB"] - ej["posB"])
            # echo(dist1, dist2)
            if dist1 < 250 and dist2 < 250 and is_reciprocal_overlapping(ei["posA"], ei["posB"], ej["posA"], ej["posB"]):
                loci_similar = True
            if dist1 < 5 and dist2 < 5:
                loci_same = True

        if not loci_similar:  # Try chrA matches chrB
            if ei["chrA"] == ej["chrB"] and ei["chrB"] == ej["chrA"]:  # Try chrA matches chrA
                dist1 = abs(ei["posA"] - ej["posB"])
                dist2 = abs(ei["posB"] - ej["posA"])
                if dist1 < 250 and dist2 < 250 and is_reciprocal_overlapping(ei["posA"], ei["posB"],
                                                                             ej["posA"], ej["posB"]):
                    loci_similar = True
                if dist1 < 5 and dist2 < 5:
                    loci_same = True

        if "contig" not in ei or "contig" not in ej:
            continue

        ci = "" if ei["contig"] is None else ei["contig"]
        ci2 = "" if ei["contig2"] is None else ei["contig2"]
        cj = "" if ej["contig"] is None else ej["contig"]
        cj2 = "" if ej["contig2"] is None else ej["contig2"]
        any_contigs_to_check = any((ei["contig"], ei["contig2"])) and any((ej["contig"], ej["contig2"]))
        # echo("loci", loci_similar, loci_same, any_contigs_to_check, ei["posA"], ei["posB"], ej["posA"], ej["posB"])

        if loci_similar or loci_same:
            if any_contigs_to_check and loci_same:  # Try and match contigs
                G.add_edge(i_id, j_id, loci_same=loci_same)
                # echo(i_id, j_id)
                # echo("loci same", ei["posA"], ei["posB"], ej["posA"], ej["posB"])

            # No contigs to match, merge anyway
            elif not any_contigs_to_check and loci_same:
                G.add_edge(i_id, j_id, loci_same=loci_same)
                # echo("2", i_id, j_id)
            # Only merge loci if they are not both within --include regions. chrA:posA only needs checking
            elif not (io_funcs.intersecter_str_chrom(tree, ei["chrA"], ei["posA"], ei["posA"] + 1) and
                      io_funcs.intersecter_str_chrom(tree, ei["chrB"], ei["posB"], ei["posB"] + 1)):
                G.add_edge(i_id, j_id, loci_same=loci_same)
                # echo("3", i_id, j_id)

            # echo(io_funcs.intersecter_str_chrom(tree, ei["chrA"], ei["posA"], ei["posA"] + 1),
            #           io_funcs.intersecter_str_chrom(tree, ei["chrB"], ei["posB"], ei["posB"] + 1))
        # echo("len g nodes", len(G.nodes()))
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
                 rel_diffs=False, diffs=15, same_sample=True, debug=False, min_size=0):
    """Try and merge similar events, use overlap of both breaks points
    """

    max_dist = max_dist / 2
    if len(potential) <= 1:
        return potential

    # Cluster events on graph
    G = nx.Graph()
    G = enumerate_events(G, potential, max_dist, try_rev, tree, rel_diffs, diffs, same_sample,
                         debug=debug)

    found = []
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item["event_id"]):
            found.append(item)

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
            svlen = w0["svlen"]
            svt = w0["svtype"]
            for k in range(1, len(best)):

                item = best[k]

                # Sum these
                for t in ("pe", "supp", "sc", "su", "NP", "block_edge", "plus", "minus", "spanning"):
                    if t in item:
                        best[0][t] += item[t]

                if item["maxASsupp"] > w0["maxASsupp"]:
                    best[0]["maxASsupp"] = item["maxASsupp"]

                if item["svtype"] == "DEL" and svt == "DEL" and (item["svlen"] * 0.6 < svlen < item["svlen"] or min_size > svlen < item["svlen"]):  # increase svlen size
                    svlen = item["svlen"]
                    best[0]["svlen"] = svlen

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


def elide_insertions(events, max_dist=150):
    return events
    # If an imprecise insertion overlaps the edge of a deletion; remove the insertion

    interval_table = defaultdict(lambda: graph.Table())
    # to_check = []
    # Make a nested containment list for lookups
    for idx in range(len(events)):
        ei = events[idx]
        if ei["svtype"] == "INS" and ei["svlen_precise"] == 0:
            interval_table[ei["chrA"]].add(ei["posA"] - ei["cipos95A"] - max_dist, ei["posA"] + ei["cipos95A"] + max_dist, idx)
            if ei["posA"] != ei["posB"]:
                dist1 = abs(ei["posB"] - ei["posA"])
                ci_a = max(ei["cipos95A"], ei["cipos95A"])
                if dist1 + ci_a > max_dist:
                    interval_table[ei["chrB"]].add(ei["posB"] - ei["cipos95B"] - max_dist, ei["posB"] + ei["cipos95B"] + max_dist, idx)
        # else:
        #     to_check.append(idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}
    if not nc:
        return events
    bad_ids = set([])

    for idx in events:
        ei = events[idx]
        ols = []

        if ei["chrA"] in nc:
            ols += list(nc[ei["chrA"]].find_overlap(ei["posA"], ei["posA"] + 1))
        if ei["chrB"] in nc:
            ols += list(nc[ei["chrB"]].find_overlap(ei["posB"], ei["posB"] + 1))

        for target in ols:
            other = events[target[2]]
            bad_ids.add(other["event_id"])

    return [i for i in events if i["event_id"] not in bad_ids]


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


cdef bint poly_kmer(k):
    cdef int len_k = len(k)
    if len_k == 1:
        return False
    sk = set(k)
    if len(sk) == 1:
        return True
    elif len_k == 4:
        if k[:2] == k[2:]:
            return True
    elif len_k == 6:
        if k[:3] == k[3:]:
            return True
    return False


cdef dict search_ssr_kc(ori):

    seq = ori.upper()
    cdef const unsigned char[:] string_view = bytes(seq.encode("ascii"))  # use for indexing

    cdef int str_len = len(seq)
    cdef int rep_len = min(7, str_len)
    cdef int i = 0
    cdef int t, start, count, mm, good_i, successive_bad, size, finish

    cdef int n_ref_repeat_bases = 0
    cdef int n_expansion = 0
    cdef int stride = 0
    expansion_seq = ""

    while i < str_len:

        if seq[i] == b"N":
            i += 1
            continue

        for t in range(1, rep_len):

            start = i
            if start + t >= str_len:
                break

            starting_base = seq[i]
            starting_kmer = seq[start:start + t]

            if poly_kmer(starting_kmer):
                continue

            count = 1
            mm = 0
            good_i = 0
            successive_bad = 0
            finish = 0
            while start + t < str_len and (starting_base == seq[start+t] or mm < 2):
                start += t
                if seq[start: start + t] != starting_kmer:
                    successive_bad += 1
                    mm += 1
                    if mm > 3 or successive_bad > 1 or count < 2:
                        finish = good_i
                        break
                else:
                    good_i = start
                    finish = good_i
                    successive_bad = 0
                    count += 1

            if count >= 3 and (finish - i) + t > 10:

                # check for lowercase to uppercase transition; determines repeat expansion length
                lowercase = []
                uppercase = []
                expansion_index = -1
                for j in range(i, finish, t):  # only check first base of repeat motif
                    if ori[j].islower():
                        if not lowercase:  # finds first block of repeats
                            lowercase.append([j, j + t])
                            if uppercase and abs(uppercase[-1][1] - j) < 3:
                                expansion_index = 0
                        elif lowercase[-1][1] == j:
                            lowercase[-1][1] += t
                    else:  # finds all reference blocks
                        if not uppercase:  # finds first block of repeats
                            uppercase.append([j, j + t])
                            if lowercase and abs(lowercase[-1][1] - j) < 3:
                                expansion_index = len(lowercase) - 1

                        elif uppercase[-1][1] == j:
                            uppercase[-1][1] += t

                        else:
                            uppercase.append([j, j + t])

                if expansion_index != -1:
                    e = lowercase[expansion_index]
                    size = e[1] - e[0]
                    if size >= 10:
                        n_expansion = size
                        expansion_seq = ori[e[0]:e[1]]
                        stride = t

                for begin, end in uppercase:
                    n_ref_repeat_bases += end - begin

                i = finish + t
        i += 1

    return {"n_expansion": n_expansion, "stride": stride, "exp_seq": expansion_seq, "ref_poly_bases": n_ref_repeat_bases}


cdef find_repeat_expansions(events, insert_stdev):
    for e in events:
        res = {"n_expansion": 0, "stride": 0, "exp_seq": "", "ref_poly_bases": 0}
        if "contig" in e and e["contig"] is not None:
            res.update(search_ssr_kc(e["contig"]))
        if "contig2" in e and e["contig2"] is not None:
            r = search_ssr_kc(e["contig2"])
            if res["n_expansion"] < r["n_expansion"]:
                res["n_expansion"] = r["n_expansion"]
                res["stride"] = r["stride"]
                res["exp_seq"] = r["exp_seq"]
            res["ref_poly_bases"] += r["ref_poly_bases"]
        e.update(res)

        # if e["svtype"] == "INS" and e["n_expansion"] > 0 and e['svlen_precise'] == 0:
        #     if abs(e["svlen"]) - insert_stdev < 0:
        #         e["svlen"] = e["n_expansion"]

    return events


def component_job(infile, component, regions, event_id, clip_length, insert_med, insert_stdev, min_supp,
                  merge_dist, regions_only, extended_tags, assemble_contigs, rel_diffs, diffs, min_size):

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
                                    try_rev=False, pick_best=True, rel_diffs=rel_diffs, diffs=diffs, min_size=min_size)

    return potential_events, event_id


def pipe1(args, infile, kind, regions, ibam, ref_genome):

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
        read_len = genome_scanner.approx_read_length
        max_dist = int(insert_median + (insert_stdev * 5))  # 5
        max_clust_dist = 1 * (int(insert_median + (5 * insert_stdev)))
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
    lower_bound_support = min_support #- 1 if min_support - 1 > 1 else 1
    echo("lower bound suport", lower_bound_support)

    read_buffer = genome_scanner.read_buffer

    t5 = time.time()
    G, node_to_name = graph.construct_graph(genome_scanner,
                                            infile,
                                            max_dist=max_dist,
                                            clustering_dist=max_clust_dist,
                                            minimizer_dist=150,
                                            minimizer_support_thresh=args["z_depth"],
                                            minimizer_breadth=args["z_breadth"],
                                            k=12,
                                            m=6,
                                            clip_l=args["clip_length"],
                                            min_sv_size=args["min_size"],
                                            min_support=lower_bound_support, #min_support,
                                            procs=args["procs"],
                                            mapq_thresh=args["mq"],
                                            debug=None,
                                            paired_end=paired_end,
                                            read_length=read_len,
                                            )
    echo("graph time", time.time() - t5)
    click.echo("graph time, mem={} Mb, time={} h:m:s".format(
        int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
        time.time() - t5), err=True)

    t0 = time.time()
    cmp, cc = graph.get_block_components(G, node_to_name, infile, read_buffer)

    extended_tags = genome_scanner.extended_tags

    if paired_end:
        rel_diffs = False
        diffs = 15
    else:
        rel_diffs = True
        diffs = 0.15
    echo("len components", len(cc))
    for start_i, end_i in cc:

        if end_i - start_i >= lower_bound_support:
            component = list(cmp[start_i: end_i])

            res = graph.proc_component(node_to_name, component, read_buffer, infile, G, lower_bound_support)
            if res:
                # Res is a dict
                # {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}
                potential_events, event_id = component_job(infile, res, regions, event_id, args["clip_length"],
                                                 insert_median,
                                                 insert_stdev,
                                                 lower_bound_support, #args["min_support"],
                                                 args["merge_dist"],
                                                 regions_only,
                                                 extended_tags,
                                                 assemble_contigs,
                                                 rel_diffs=rel_diffs, diffs=diffs, min_size=args["min_size"])
                if potential_events:
                    block_edge_events += potential_events
                    # echo(potential_events)

    echo("components", time.time() - t0)
    t3 = time.time()
    preliminaries = []

    tremap = time.time()
    block_edge_events = re_map.remap_soft_clips(block_edge_events, ref_genome, args["min_size"])
    echo("remapping", time.time() - tremap)
    # for item in block_edge_events:
    #     echo(item)
    # Merge across calls
    if args["merge_within"] == "True":
        tmr = time.time()
        merged = merge_events(block_edge_events, args["merge_dist"], regions, try_rev=False, pick_best=False,
                                         debug=True, min_size=args["min_size"])
        echo("merging all", time.time() - tmr)
    else:
        merged = block_edge_events

    merged = elide_insertions(merged)

    # Filter for absolute support and size here
    merged = [event for event in merged if event["svlen"] >= args["min_size"] and event["su"] >= args["min_support"]]
    # for e in merged:
    #     if e["event_id"] == 6:
    #         echo(e)
    if merged:
        for event in merged:
            # Collect coverage information
            event_dict = call_component.get_raw_coverage_information(event, regions, genome_scanner.depth_d, infile) # regions_depth)

            if event_dict:
                preliminaries.append(event_dict)

    preliminaries = assembler.contig_info(preliminaries)  # GC info, repetitiveness
    preliminaries = sample_level_density(preliminaries, regions)
    preliminaries = find_repeat_expansions(preliminaries, insert_stdev)

    if len(preliminaries) == 0:
        click.echo("No events found", err=True)
        quit()
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

    ref_genome = pysam.FastaFile(args["reference"])

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
    events, extended_tags = pipe1(args, infile, kind, regions, ibam, ref_genome)
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
    ref_genome.close()

    echo(time.time() - t0)
    click.echo("dysgu call {} complete, n={}, time={} h:m:s".format(
               args["sv_aligns"],
               len(df),
               str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)
