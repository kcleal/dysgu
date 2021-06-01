# cython: language_level=3

from __future__ import absolute_import
import datetime
import os
import time
import logging
import numpy as np
import random
from collections import defaultdict, Counter
import networkx as nx
import pysam
import sys
import resource
import pandas as pd
from dysgu import coverage, graph, call_component, assembler, io_funcs, re_map, post_call_metrics
from dysgu.map_set_utils cimport is_reciprocal_overlapping, Py_CoverageTrack, EventResult
from dysgu.map_set_utils import timeit, echo
import itertools


import itertools
import multiprocessing
from scipy import stats
from libcpp.vector cimport vector


ctypedef EventResult EventResult_t


def filter_potential(input_events, tree, regions_only):
    potential = []
    cdef EventResult_t i
    for i in input_events:

        # if "posB" not in i:  # Skip events for which no posB was identified
        #     continue
        if i.sqc is None:
            i.sqc = -1
        if i.svtype == "INS" and i.svlen_precise == 0 and not i.contig and not i.contig2:
            continue
        # if "contig" not in i or i["contig"] == "":
        #     i["contig"] = None

        # if i.bnd is None:
        #     i.bnd = 0
        # if "remap_score" not in i or i["remap_score"] is None:
        #     i["remap_score"] = 0
        # if "scw" not in i or i["scw"] is None:
        #     i["scw"] = 0
        # Remove events for which both ends are in --regions but no contig was found
        posA_intersects = io_funcs.intersecter_str_chrom(tree, i.chrA, i.posA, i.posA + 1)
        posB_intersects = io_funcs.intersecter_str_chrom(tree, i.chrB, i.posB, i.posB + 1)

        if (posA_intersects and posB_intersects) and (i.contig is None or i.contig2 is None):
            continue
        # Remove events for which neither end is in --regions (if --regions provided)
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
        interval_table[ei.chrA].add(ei.posA - ei.cipos95A - max_dist, ei.posA + ei.cipos95A + max_dist, idx)
        # Add another interval for large events, or translocations
        if ei.chrA != ei.chrB or abs(ei.posB - ei.posA) > 5:
            if ei.chrA != ei.chrB:
                interval_table[ei.chrB].add(ei.posB - ei.cipos95B - max_dist, ei.posB + ei.cipos95B + max_dist, idx)
            else:
                dist1 = abs(ei.posB - ei.posA)
                ci_a = max(ei.cipos95A, ei.cipos95A)
                if dist1 + ci_a > max_dist:
                    interval_table[ei.chrB].add(ei.posB - ei.cipos95B - max_dist, ei.posB + ei.cipos95B + max_dist, idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}

    for idx in range(len(potential)):
        ei = potential[idx]
        # Overlap right hand side
        ols = list(nc[ei.chrB].find_overlap(ei.posB, ei.posB + 1))
        for target in ols:
            jdx = target[2]
            ej = potential[jdx]
            yield ei, ej, idx, jdx


def compare_all(potential):
    # Quadratic run time but fast for a hand full of events
    for idx, jdx in itertools.product(range(len(potential)), range(len(potential))):
        yield potential[idx], potential[jdx], idx, jdx


cdef span_similarity(EventResult_t ei, EventResult_t ej):
    if ei.svtype != "INS":
        span1 = abs(ei.posB - ei.posA)
        span2 = abs(ej.posB - ej.posA)
        max_span = max(span1, span2)
        if max_span > 0:
            span_distance = abs(span1 - span2) / max_span
            return span_distance < 0.1
    return False


cdef break_distances(i_a, i_b, j_a, j_b, i_a_precise, i_b_precise, j_a_precise, j_b_precise, intra=True):

    if intra:
        if i_b < i_a:
            temp = i_a
            i_a = i_b
            i_b = temp
            temp = i_a_precise
            i_a_precise = i_b_precise

        if j_b < j_a:
            temp = j_a
            j_a = j_b
            j_b = temp
            temp = j_a_precise
            j_a_precise = j_b_precise

    dist_a = abs(i_a - j_a)
    dist_b = abs(i_b - j_b)

    precise_thresh = 250
    imprecise_thresh = 650
    same_thresh = 5

    dist_a_same = dist_a < same_thresh
    if not i_a_precise or not j_a_precise:
        dist_a_close = dist_a < imprecise_thresh
    else:
        dist_a_close = dist_a < precise_thresh

    dist_b_same = dist_b < same_thresh
    if not i_b_precise or not j_b_precise:
        dist_b_close = dist_b < imprecise_thresh
    else:
        dist_b_close = dist_b < precise_thresh

    return dist_a_close and dist_b_close, dist_a_same and dist_b_same


def enumerate_events(G, potential, max_dist, try_rev, tree, paired_end=False, rel_diffs=False, diffs=15,
                     same_sample=True, aggressive_ins_merge=False, debug=False):

    if len(potential) < 3:
        event_iter = compare_all(potential)  # N^2 compare all events to all others
    else:
        event_iter = compare_subset(potential, max_dist)  # Use NCLS, generate overlap tree and perform intersections

    seen = set([])
    pad = 100
    # cdef EventResult_t ei, ej
    for ei, ej, idx, jdx in event_iter:

        i_id = ei.event_id
        j_id = ej.event_id
        if not same_sample: # and "sample" in ei:
            if ei.sample == ej.sample:
                continue
        else:
            if i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen:
                continue

        seen.add((i_id, j_id))

        fail = False
        if ei.svtype != ej.svtype:
            fail = True
        if fail:
            continue

        # Check if events point to the same loci
        both_in_include = io_funcs.intersecter_str_chrom(tree, ei.chrA, ei.posA, ei.posA + 1) and \
                      io_funcs.intersecter_str_chrom(tree, ei.chrB, ei.posB, ei.posB + 1)
        loci_similar = False
        loci_same = False
        intra = ei.chrA == ei.chrB == ej.chrA == ej.chrB

        if intra:
            loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posA, ej.posB, ei.preciseA, ei.preciseB, ej.preciseA, ej.preciseB)

        elif ei.chrA == ej.chrA and ei.chrB == ej.chrB:  # Try chrA matches chrA
            loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posA, ej.posB, ei.preciseA, ei.preciseB, ej.preciseA, ej.preciseB, intra=False)

        elif ei.chrA == ej.chrB and ei.chrB == ej.chrA:
            loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posB, ej.posA, ei.preciseA, ei.preciseB, ej.preciseB, ej.preciseA, intra=False)

        one_is_imprecise = False
        if ((not ei.preciseA or not ei.preciseB) and ei.svlen_precise) or ((not ej.preciseA or not ej.preciseB) and ej.svlen_precise):
            one_is_imprecise = True  # also not remapped

        if not loci_similar:
            continue

        if loci_same:
            if not both_in_include:
                G.add_edge(i_id, j_id, loci_same=loci_same)
                continue

        ci = ei.contig
        ci2 = ei.contig2

        cj = ej.contig
        cj2 = ej.contig2

        any_contigs_to_check = any((ci, ci2)) and any((cj, cj2))
        recpi_overlap = is_reciprocal_overlapping(ei.posA, ei.posB, ej.posA, ej.posB)

        # If long reads only rely on reciprocal overlap, seems to work better
        if paired_end:
            spd = span_similarity(ei, ej)
        else:
            spd = False

        m = False
        ml = max(ei.svlen, ej.svlen)
        if ei.svtype == "INS":
            if aggressive_ins_merge:
                m = True
            else:
                # if both have been remapped, make sure size is similar
                if ei.spanning > 0 or ej.spanning > 0:
                    if paired_end:
                        m = True
                    elif ml > 0 and min(ei.svlen, ej.svlen) / ml > 0.8:
                        m = True
                elif ei.remap_score > 0 and ej.remap_score > 0 and ml > 0 and min(ei.svlen, ej.svlen) / ml > 0.8:
                    m = True
                elif ei.remap_score == 0 or ej.remap_score == 0:  # merge if one break was not remapped
                    m = True

        elif ei.svtype != "INS" and (recpi_overlap or spd) and not both_in_include:
            m = True

        if not m:
            continue

        # Loci are similar, check contig match or reciprocal overlap
        if not any_contigs_to_check:
            if ml > 0:
                l_ratio = min(ei.svlen, ej.svlen) / ml
                if l_ratio > 0.5 or (one_is_imprecise and l_ratio > 0.3):
                    G.add_edge(i_id, j_id, loci_same=loci_same)

        else:
            if ci and cj: v = ci, cj
            elif ci and cj2: v = ci, cj2
            elif ci2 and cj2: v = ci2, cj2
            elif ci2 and cj: v = ci2, cj
            else: continue
            # if not remapped and nearby insertions with opposing soft-clips --> merge
            # also long-reads will normally have remap_score == 0
            if ei.remap_score == 0 or ej.remap_score == 0:
                if (v[0][0].islower() and v[1][-1].islower()) or (v[0][-1].islower() and v[1][0].islower()):
                    G.add_edge(i_id, j_id, loci_same=loci_same)
                    continue

            if assembler.check_contig_match(v[0], v[1], return_int=True):
                G.add_edge(i_id, j_id, loci_same=True)

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
    # if "type" in c and c["type"] != "1":  # assume long read is != "1", either pacbio or nanopore
    if c.type != "pe" and c.type != "":
        return 100 + c.su
    # if "su" in c:
    return c.su
    # over weight spanning
    # return c["pe"] + c["supp"] + (c["spanning"] * 10) + c["sc"]


def merge_events(potential, max_dist, tree, paired_end=False, try_rev=False, pick_best=False, add_partners=False,
                 rel_diffs=False, diffs=15, same_sample=True, debug=False, min_size=0, aggressive_ins_merge=False,
                 skip_imprecise=False):
    """Try and merge similar events, use overlap of both breaks points
    """

    max_dist = max_dist / 2
    if len(potential) <= 1:
        return potential

    # Cluster events on graph
    G = nx.Graph()
    G = enumerate_events(G, potential, max_dist, try_rev, tree, paired_end, rel_diffs, diffs, same_sample,
                         aggressive_ins_merge=aggressive_ins_merge,
                         debug=debug)

    found = []
    # cdef EventResult_t item, w0
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item.event_id):
            found.append(item)

    # Try and merge SVs with identical breaks, then merge ones with less accurate breaks - this helps prevent
    # over merging SVs that are close together
    components = cut_components(G)
    node_to_event = {i.event_id: i for i in potential}
    cdef int k
    # Only keep edges with loci_same==False if removing the edge leads to an isolated node
    for grp in components:

        c = [node_to_event[n] for n in grp]
        best = sorted(c, key=srt_func, reverse=True)
        w0 = best[0]
        if not pick_best:

            # Weighting for base result
            weight = w0.pe + w0.supp + w0.spanning
            spanned = bool(w0.spanning)
            if w0.svlen_precise >= 0:
                remapped = bool(not w0.svlen_precise and w0.remap_score > 0)
            else:
                remapped = False
            add_contig_a = bool(not w0.contig)
            new_a = ""
            add_contig_b = bool(not w0.contig2)
            new_b = ""
            svt = w0.svtype
            sqc = w0.sqc #w0["sqc"] if "sqc" in w0 else -1
            for k in range(1, len(best)):

                item = best[k]
                # Sum these
                w0.pe += item.pe
                w0.supp += item.supp
                w0.sc += item.sc
                w0.su += item.su
                w0.NP += item.NP
                w0.block_edge += item.block_edge
                w0.plus += item.plus
                w0.minus += item.minus
                w0.spanning += item.spanning
                w0.n_small_tlen += item.n_small_tlen

                if item.maxASsupp > w0.maxASsupp:
                    w0.maxASsupp = item.maxASsupp

                if item.svtype == "DEL" and svt == "DEL":
                    if not spanned:

                        if item.spanning:
                            w0.svlen = item.svlen

                        elif min_size > w0.svlen < item.svlen:  # increase svlen size
                            w0.svlen = item.svlen

                elif item.svtype == "INS" and svt == "INS":
                    if not spanned:
                        if item.spanning:
                            w0.svlen = item.svlen

                        elif item.svlen * 0.6 < w0.svlen < item.svlen or min_size > w0.svlen < item.svlen:  # increase svlen size
                            w0.svlen = item.svlen

                if item.sqc != -1 and sqc == 1:
                    w0.sqc = sqc
                    sqc = item.sqc

                # Average weight these
                wt = item.pe + item.supp + item.spanning
                denom = weight + wt

                def norm_vals(a, wt_a, b, wt_b, denom):
                    if denom > 0 and a is not None and b is not None:
                        return ((a * wt_a) + (b * wt_b)) / denom
                    return a

                if denom > 0:
                    w0.DN = norm_vals(w0.DN, weight, item.DN, wt, denom)
                    w0.MAPQsupp = norm_vals(w0.MAPQsupp, weight, item.MAPQsupp, wt, denom)
                    w0.MAPQpri = norm_vals(w0.MAPQpri, weight, item.MAPQpri, wt, denom)
                    w0.DApri = norm_vals(w0.DApri, weight, item.DApri, wt, denom)
                    w0.DAsupp = norm_vals(w0.DAsupp, weight, item.DAsupp, wt, denom)
                    w0.DP = norm_vals(w0.DP, weight, item.DP, wt, denom)
                    w0.NMpri = norm_vals(w0.NMpri, weight, item.NMpri, wt, denom)
                    w0.NMsupp = norm_vals(w0.NMsupp, weight, item.NMsupp, wt, denom)

                    # w0.DN = ((w0.DN * weight) + (item.DN * wt)) / denom
                    # w0.MAPQsupp = ((w0.MAPQsupp * weight) + (item.MAPQsupp * wt)) / denom
                    # w0.MAPQpri = ((w0.MAPQpri * weight) + (item.MAPQpri * wt)) / denom
                    # w0.DApri = ((w0.DApri * weight) + (item.DApri * wt)) / denom
                    # w0.DAsupp = ((w0.DAsupp * weight) + (item.DAsupp * wt)) / denom
                    # w0.DP = ((w0.DP * weight) + (item.DP * wt)) / denom
                    # w0.NMpri = ((w0.NMpri * weight) + (item.NMpri * wt)) / denom
                    # w0.NMsupp = ((w0.NMsupp * weight) + (item.NMsupp * wt)) / denom

                # reset to combined weight?
                # weight = w0.pe + w0.supp + w0.spanning

                if add_contig_a and item.contig and len(item.contig) > len(new_a):
                    new_a = item.contig
                if add_contig_b and item.contig2 and len(item.contig2) > len(new_b):
                    new_b = item.contig2

            if add_contig_a and new_a:
                w0.contig = new_a
            if add_contig_b and new_b:
                w0.contig2 = new_b

        if add_partners:
            # echo([i.event_id for i in best[1:]])
            w0.partners = [i.event_id for i in best[1:]]

        found.append(w0)

    return found


def sample_level_density(potential, regions, max_dist=50):

    interval_table = defaultdict(lambda: graph.Table())

    # Make a nested containment list for faster lookups
    cdef EventResult_t ei
    for idx in range(len(potential)):
        ei = potential[idx]

        # Only find density for non-region calls, otherwise too dense to be meaningful
        if not io_funcs.intersecter_str_chrom(regions, ei.chrA, ei.posA, ei.posA + 1):
            interval_table[ei.chrA].add(ei.posA - max_dist, ei.posA + max_dist, idx)

        if not io_funcs.intersecter_str_chrom(regions, ei.chrB, ei.posB, ei.posB + 1):
            interval_table[ei.chrB].add(ei.posB - max_dist, ei.posB + max_dist, idx)

    nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}

    for idx in range(len(potential)):
        ei = potential[idx]

        # Overlap right hand side # list(tree[chrom].find_overlap(start, end))
        neighbors = 0.
        count = 0.
        if ei.chrA == ei.chrB and abs(ei.posB - ei.posA) < 2:
            expected = 2
        else:
            expected = 1
        if not io_funcs.intersecter_str_chrom(regions, ei.chrA, ei.posA, ei.posA + 1):
            neighbors += len(list(nc[ei.chrA].find_overlap(ei.posA, ei.posA + 1))) - expected
            count += 1

        if not io_funcs.intersecter_str_chrom(regions, ei.chrB, ei.posB, ei.posB + 1):
            neighbors += len(list(nc[ei.chrB].find_overlap(ei.posB, ei.posB + 1))) - expected
            count += 1

        # Merge overlapping intervals
        neighbors_10kb = 0.
        count_10kb = 0
        large_itv = coverage.merge_intervals(((ei.chrA, ei.posA, ei.posA + 1),
                                              (ei.chrB, ei.posB, ei.posB + 1)), pad=10000)

        for c, s, e in large_itv:
            if not io_funcs.intersecter_str_chrom(regions, c, s, e):
                neighbors_10kb += len(list(nc[c].find_overlap(s, e))) - len(large_itv)
                count_10kb += 1

        if neighbors < 0:
            neighbors = 0
        if count > 0:
            ei.neigh = neighbors / count
        else:
            ei.neigh = 0
        if count_10kb > 0:
            ei.neigh10kb = neighbors_10kb / count_10kb
        else:
            ei.neigh10kb = 0

    return potential


cdef bint same_k(int start1, int start2, int n, const unsigned char[:] seq):
    cdef int j
    for j in range(n):
        if seq[start1 + j] != seq[start2 + j]:
            return False
    return True


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

    cdef unsigned char starting_base
    while i < str_len:

        if seq[i] == b"N":
            i += 1
            continue

        for t in range(1, rep_len):

            start = i
            if start + t >= str_len:
                break

            starting_base = string_view[i]
            starting_kmer_idx = start

            # if seq[start:start+t] in polymers:
            #     continue

            count = 1
            mm = 0
            good_i = 0
            successive_bad = 0
            finish = 0
            while start + t < str_len and (starting_base == string_view[start+t] or mm < 2):
                start += t
                if start + t + 1 > str_len:
                    break

                if not same_k(start, starting_kmer_idx, t, string_view):
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


# @timeit
def find_repeat_expansions(events, insert_stdev):

    cdef EventResult_t e
    for e in events:
        e.n_expansion = 0
        e.stride = 0
        e.exp_seq = ""
        e.ref_poly_bases = 0
        if e.contig:
            r = search_ssr_kc(e.contig)
            e.n_expansion = r["n_expansion"]
            e.stride = r["stride"]
            e.exp_seq = r["exp_seq"]
            e.ref_poly_bases += r["ref_poly_bases"]

        if e.contig2:
            r = search_ssr_kc(e.contig2)
            if e.n_expansion < r["n_expansion"]:
                e.n_expansion = r["n_expansion"]
                e.stride = r["stride"]
                e.exp_seq = r["exp_seq"]
                e.ref_poly_bases += r["ref_poly_bases"]

    return events


def component_job(infile, component, regions, event_id, clip_length, insert_med, insert_stdev, insert_ppf, min_supp, lower_bound_support,
                  merge_dist, regions_only, extended_tags, assemble_contigs, rel_diffs, diffs, min_size,
                  ):

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
                                                      extended_tags,
                                                      assemble_contigs,
                                                      ):

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


class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, bam_path, open_mode, regions, clip_length, insert_median, insert_stdev, insert_ppf, support, lower_bound_support, merge_dist,
                 regions_only, extended_tags, assemble_contigs, rel_diffs, diffs, min_size):
        self.infile = pysam.AlignmentFile(bam_path, open_mode)
        self.regions = regions
        self.clip_length = clip_length
        self.insert_median = insert_median
        self.insert_stdev = insert_stdev
        self.insert_ppf = insert_ppf
        self.support = support
        self.lower_bound_support = lower_bound_support
        self.merge_dist = merge_dist
        self.regions_only = regions_only
        self.extended_tags = extended_tags
        self.assemble_contigs = assemble_contigs
        self.rel_diffs = rel_diffs
        self.diffs = diffs
        self.min_size = min_size
        # these must be at the bottom:
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.infile.close()
                self.task_queue.task_done()
                break
            event_id, job = next_task()
            events = []
            potential_events, event_id = component_job(self.infile,
                                                       job,
                                                       self.regions,
                                                       event_id,
                                                       self.clip_length,
                                                       self.insert_median,
                                                       self.insert_stdev,
                                                       self.insert_ppf,
                                                       self.support,
                                                       self.lower_bound_support,
                                                       self.merge_dist,
                                                       self.regions_only,
                                                       self.extended_tags,
                                                       self.assemble_contigs,
                                                       self.rel_diffs,
                                                       self.diffs,
                                                       self.min_size)

            if potential_events:
                events += potential_events
            self.result_queue.put(events)
            self.task_queue.task_done()


class Task:
    def __init__(self, job_id, job):
        self.job_id = job_id
        self.job = job
    def __call__(self):
        return self.job_id, self.job
    def __str__(self):
        return f'{self.job_id} * {len(self.job)}'


def pipe1(args, infile, kind, regions, ibam, ref_genome, open_mode):

    regions_only = False if args["regions_only"] == "False" else True
    paired_end = int(args["paired"] == "True")
    assemble_contigs = int(args["contigs"] == "True")
    temp_dir = args["working_directory"]
    if not args["ibam"]:
        # Make a new coverage track if one hasn't been created yet
        coverage_tracker = Py_CoverageTrack(temp_dir, infile, args["max_cov"])
    else:
        coverage_tracker = None

    genome_scanner = coverage.GenomeScanner(infile, args["max_cov"], args["regions"], 1,
                                            args["buffer_size"], regions_only,
                                            kind == "stdin",
                                            clip_length=args["clip_length"],
                                            min_within_size=args["min_size"],
                                            coverage_tracker=coverage_tracker)

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
        logging.info(f"Max clustering dist {max_clust_dist}")
    else:
        if args["mode"] == "pacbio":
            max_dist, max_clust_dist = 35, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 50
        elif args["mode"] == "nanopore":
            max_dist, max_clust_dist = 100, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 150

    event_id = 0
    block_edge_events = []
    min_support = args["min_support"]
    logging.info("Minimum support {}".format(min_support))
    if args["pl"] == "pe":  # reads with internal SVs can be detected at lower support
        lower_bound_support = min_support - 1 if min_support - 1 > 1 else 1
    else:
        lower_bound_support = min_support
    clip_length = args["clip_length"]
    merge_dist = args["merge_dist"]
    min_size = args["min_size"]

    read_buffer = genome_scanner.read_buffer

    t5 = time.time()
    G, node_to_name, bad_clip_counter = graph.construct_graph(genome_scanner,
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
                                            procs=1, #args["procs"],
                                            mapq_thresh=args["mq"],
                                            debug=None,
                                            paired_end=paired_end,
                                            read_length=read_len,
                                            contigs=args["contigs"],
                                            norm_thresh=args["dist_norm"],
                                            spd_thresh=args["spd"],
                                            mm_only=args["regions_mm_only"] == "True")

    logging.info("Graph constructed")
    # logging.info("Graph time, mem={} Mb, time={} h:m:s".format(
    #     int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
    #     time.time() - t5))

    t0 = time.time()
    cdef vector[int] cmp = G.connectedComponents()  # Flat vector, components are separated by -1
    # logging.info("connected components done")
    if insert_median != -1:
        insert_ppf = stats.norm.ppf(0.05, loc=insert_median, scale=insert_stdev)
        if insert_ppf < 0:
            insert_ppf = 0
    else:
        insert_ppf = -1
    extended_tags = genome_scanner.extended_tags

    if paired_end:
        rel_diffs = False
        diffs = 15
    else:
        rel_diffs = True
        diffs = 0.15

    if args["procs"] > 1:
        tasks = multiprocessing.JoinableQueue(maxsize=100)
        results_q = multiprocessing.Queue(maxsize=10000)
        consumers = [ Consumer(tasks, results_q, args["sv_aligns"], open_mode, regions, clip_length, insert_median,
                               insert_stdev, insert_ppf, min_support, merge_dist, regions_only, extended_tags, assemble_contigs,
                               rel_diffs, diffs, min_size)
            for i in range(args["procs"])]

        for w in consumers:
            w.start()

    else:
        tasks = None
        results = None
        consumers = None

    num_jobs = 0
    completed = 0
    components_seen = 0
    cdef int last_i = 0
    cdef int ci, cmp_idx
    # longest = 0
    for item_idx, item in enumerate(cmp):
        if item == -1:
            components_seen += 1
            start_i = last_i
            end_i = item_idx
            last_i = item_idx + 1
            # if end_i - start_i > longest:
            #     longest = end_i - start_i

            # todo dont copy this. use memory view slice? sending vector[int]& seems to result in copying (really slow)
            component = np.zeros(end_i - start_i)
            ci = 0
            for cmp_idx in range(start_i, end_i):
                component[ci] = cmp[cmp_idx]
                ci += 1

            res = graph.proc_component(node_to_name, component, read_buffer, infile, G, lower_bound_support, args["procs"], paired_end)

            if res:
                event_id += 1
                # Res is a dict
                # {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}
                if args["procs"] == 1:

                    potential_events, event_id = component_job(infile, res, regions, event_id, clip_length,
                                                               insert_median,
                                                               insert_stdev,
                                                               insert_ppf,
                                                               min_support,
                                                               lower_bound_support,
                                                               merge_dist,
                                                               regions_only,
                                                               extended_tags,
                                                               assemble_contigs,
                                                               rel_diffs=rel_diffs, diffs=diffs, min_size=min_size)

                    if potential_events:
                        block_edge_events += potential_events

                else:

                    while not results_q.empty():
                        item = results_q.get()
                        block_edge_events += item
                        completed += 1
                        num_jobs -= 1

                    tasks.put(Task(event_id, res))
                    num_jobs += 1

    if args["procs"] > 1:

        for _ in range(args["procs"]):
            tasks.put(None)
        tasks.join()

        while num_jobs:
            item = results_q.get()
            block_edge_events += item
            completed += 1
            num_jobs -= 1

    logging.info("Number of components {}. N candidates {}".format(components_seen, len(block_edge_events)))

    del G
    del read_buffer
    cmp.clear()
    preliminaries = []
    if args["remap"] == "True" and args["contigs"] == "True":
        block_edge_events = re_map.remap_soft_clips(block_edge_events, ref_genome,
                                                    keep_unmapped=True if args["pl"] == "pe" else False,
                                                    min_support=args["min_support"])
        logging.info("Re-alignment of soft-clips done. N candidates {}".format(len(block_edge_events)))


    # Merge across calls
    if args["merge_within"] == "True":
        merged = merge_events(block_edge_events, args["merge_dist"], regions, paired_end, try_rev=False, pick_best=False,
                                         debug=True, min_size=args["min_size"])
    else:
        merged = block_edge_events
    logging.info("Number of candidate SVs merged: {}".format(len(block_edge_events) - len(merged)))
    logging.info("Number of candidate SVs after merge: {}".format(len(merged)))
    # Filter for absolute support and size here
    before = len(merged)
    if not args["keep_small"]:
        merged = [event for event in merged if (event.svlen >= args["min_size"] or event.chrA != event.chrB) and event.su >= args["min_support"]]
        logging.info("Number of candidate SVs dropped with sv-len < min-size or support < min support: {}".format(before - len(merged)))
    else:
        merged = [event for event in merged if event.su >= args["min_support"]]
        logging.info("Number of candidate SVs dropped with support < min support: {}".format(before - len(merged)))


    # Add read-type information
    cdef EventResult_t d
    for d in merged:
        d.type = args["pl"]

    merged = re_map.drop_svs_near_reference_gaps(merged, paired_end, ref_genome, args["drop_gaps"] == "True")

    coverage_analyser = post_call_metrics.CoverageAnalyser(temp_dir)
    # logging.info("Cov analyser")
    preliminaries = coverage_analyser.process_events(merged)
    # logging.info("Cov finished")
    preliminaries = coverage.get_raw_coverage_information(merged, regions, coverage_analyser, infile, args["max_cov"])  # genome_scanner.depth_d
    # logging.info("Cov info finished")
    preliminaries = assembler.contig_info(preliminaries)  # GC info, repetitiveness
    # logging.info("contig info finished")
    preliminaries = sample_level_density(preliminaries, regions)
    # logging.info("density finished")
    preliminaries = find_repeat_expansions(preliminaries, insert_stdev)
    # logging.info("repeat expansion finished")
    preliminaries = post_call_metrics.get_badclip_metric(preliminaries, bad_clip_counter, infile, regions)
    # logging.info("bcc finished")
    preliminaries = post_call_metrics.get_gt_metric(preliminaries, ibam, add_gt=args["no_gt"])
    # logging.info("gt finished")
    preliminaries = coverage_analyser.normalize_coverage_values(preliminaries)
    # logging.info("norm cov finished")

    preliminaries = post_call_metrics.compressability(preliminaries)

    preliminaries = post_call_metrics.ref_repetitiveness(preliminaries, args["mode"], ref_genome)

    n_in_grp = Counter([d.grp_id for d in preliminaries])
    for d in preliminaries:
        d.n_in_grp = n_in_grp[d.grp_id]

    if len(preliminaries) == 0:
        logging.critical("No events found")
        quit()

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

    if kind == "stdin" and args["regions"] is not None:
        raise ValueError("Cannot use stdin and include (an indexed bam is needed)")

    bam_mode = opts[kind]
    infile = pysam.AlignmentFile(args["sv_aligns"], bam_mode, threads=1)

    has_index = True
    try:
        infile.check_index()
    except ValueError:
        has_index = False
    logging.info("Input file has index {}".format(has_index))
    if not has_index and args["regions"] is not None:
        logging.info("Indexing input alignment file")
        infile.close()
        pysam.index(args["sv_aligns"])
        infile = pysam.AlignmentFile(args["sv_aligns"], bam_mode, threads=1)

    ref_genome = pysam.FastaFile(args["reference"])

    ibam = None
    if args["ibam"] is not None:
        kind = args["sv_aligns"].split(".")[-1]
        if kind == "stdin" or kind == "-" or kind not in opts:
            raise ValueError("--ibam must be a .bam/cam/sam file")
        ibam = pysam.AlignmentFile(args["ibam"], opts[kind], threads=1)

    if "RG" in infile.header:
        rg = infile.header["RG"]
        if len(rg) > 1:
            logging.warning("Warning: more than one @RG, using first sample (SM) for output: {}".format(rg[0]["SM"]))
        sample_name = rg[0]["SM"]
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
        outfile = sys.stdout
    else:
        logging.info("Writing SVs to {}".format(args["svs_out"]))
        outfile = open(args["svs_out"], "w")
    logging.info("Running pipeline")
    _debug_k = []
    regions = io_funcs.overlap_regions(args["regions"])

    # Run dysgu here:
    events, extended_tags = pipe1(args, infile, kind, regions, ibam, ref_genome, bam_mode)

    # logging.info("Extended tags: {}".format(extended_tags))

    df = pd.DataFrame.from_records([e.to_dict() for e in events])
    if not extended_tags:
        for cl in ("DN", "DP", "DApri", "DAsupp"):
            del df[cl]

    df = post_call_metrics.apply_model(df, args["pl"], args["contigs"], args["diploid"], args["paired"], args["thresholds"])

    if len(df) > 0:
        df = df.sort_values(["kind", "chrA", "posA"])
        df["sample"] = [sample_name] * len(df)
        df["id"] = range(len(df))

        df.rename(columns={"contig": "contigA", "contig2": "contigB"}, inplace=True)
        if args["out_format"] == "csv":
            cols = io_funcs.col_names(extended_tags, not args["metrics"])
            fmt = cols.pop()
            cols += fmt
            df[cols].to_csv(outfile, index=False)
        else:
            contig_header_lines = ""
            for item in infile.header["SQ"]:
                contig_header_lines += f"\n##contig=<ID={item['SN']},length={item['LN']}>"
            args["add_kind"] = "True"
            args["sample_name"] = sample_name
            io_funcs.to_vcf(df, args, {sample_name}, outfile, show_names=False, contig_names=contig_header_lines,
                            extended_tags=extended_tags)

    logging.info("dysgu call {} complete, n={}, time={} h:m:s".format(
               args["sv_aligns"],
               len(df),
               str(datetime.timedelta(seconds=int(time.time() - t0)))))
