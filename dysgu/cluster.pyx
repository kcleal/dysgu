# cython: language_level=3
from __future__ import absolute_import
import datetime
import os
import heapq
import logging
import numpy as np
import random
from collections import defaultdict, Counter
import networkx as nx
import pysam
from sys import stdout
import pandas as pd
from dysgu import coverage, graph, call_component, assembler, io_funcs, re_map, post_call
from dysgu.map_set_utils cimport is_reciprocal_overlapping, EventResult, Py_SimpleGraph
from dysgu.map_set_utils import to_dict, merge_intervals, echo
from dysgu import sites_utils
from dysgu.io_funcs import intersecter
import pickle
import gc
import itertools
import multiprocessing
from scipy import stats
from libcpp.vector cimport vector
import time

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


def compare_subset(potential, int max_dist):
    tmp_list = defaultdict(list)
    cdef int idx, jdx, dist1, ci_a
    for idx in range(len(potential)):
        ei = potential[idx]
        tmp_list[ei.chrA].append((ei.posA - ei.cipos95A - max_dist, ei.posA + ei.cipos95A + max_dist, idx))
        if ei.chrA != ei.chrB or abs(ei.posB - ei.posA) > 5:
            if ei.chrA != ei.chrB:
                tmp_list[ei.chrB].append((ei.posB - ei.cipos95B - max_dist, ei.posB + ei.cipos95B + max_dist, idx))
            else:
                dist1 = abs(ei.posB - ei.posA)
                ci_a = max(ei.cipos95A, ei.cipos95A)
                if dist1 + ci_a > max_dist:
                    tmp_list[ei.chrB].append((ei.posB - ei.cipos95B - max_dist, ei.posB + ei.cipos95B + max_dist, idx))
    nc2 = {k: io_funcs.iitree(v, add_value=True) for k, v in tmp_list.items()}
    for idx in range(len(potential)):
        ei = potential[idx]
        ols = nc2[ei.chrB].allOverlappingIntervals(ei.posB, ei.posB + 1)
        for jdx in ols:
            ej = potential[jdx]
            yield ei, ej, idx, jdx


def compare_all(potential):
    for idx, jdx in itertools.product(range(len(potential)), range(len(potential))):
        yield potential[idx], potential[jdx], idx, jdx


cdef span_position_distance(ei, ej):
    if ei.svtype != "INS":
        span1 = abs(ei.posB - ei.posA)
        span2 = abs(ej.posB - ej.posA)
        max_span = max(span1, span2)
        center1 = (ei.posB + ei.posA) / 2
        center2 = (ej.posB + ej.posA) / 2
        if max_span > 0:
            span_distance = abs(span1 - span2) / max_span
            position_distance = abs(center1 - center2)
            spd = (position_distance / 900) + span_distance
            return spd
    else:
        span1 = ei.svlen
        span2 = ej.svlen
        max_span = max(span1, span2)
        center1 = ei.posA
        center2 = ej.posA
        if max_span > 0:
            span_distance = abs(span1 - span2) / max_span
            position_distance = abs(center1 - center2)
            spd = (position_distance / 900) + span_distance
            return spd
    return 0


cdef break_distances(int i_a, int i_b, int j_a, j_b, bint i_a_precise, bint i_b_precise, bint j_a_precise,
                     bint j_b_precise, int min_svlen, bint intra=True):
    cdef int temp, dist_a, dist_b, precise_thresh, mprecise_thresh, same_thresh
    cdef bint dist_a_same, dist_a_close, dist_b_same, dist_b_close
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
    precise_thresh = min(max(350, min_svlen), 1500)
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
    if len(potential) < 50:
        event_iter = compare_all(potential)  # N^2 compare all events to all others
    else:
        event_iter = compare_subset(potential, max_dist)  # Use iitree, generate overlap tree and perform intersections
    seen = set([])
    pad = 100
    disjoint_nodes = set([])  # if a component has more than one disjoint nodes it needs to be broken apart
    for ei, ej, idx, jdx in event_iter:
        i_id = ei.event_id
        j_id = ej.event_id
        if not same_sample:
            if ei.sample == ej.sample:
                continue
        if i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen:
            continue
        seen.add((i_id, j_id))
        fail = False
        if ei.svtype == "INS" or ei.svtype == "DUP":
            if not (ej.svtype == "INS" or ej.svtype == "DUP"):
                continue
        elif ei.svtype != ej.svtype:
            continue

        # Check if events point to the same loci
        both_in_include = intersecter(tree, ei.chrA, ei.posA, ei.posA + 1) and intersecter(tree, ei.chrB, ei.posB, ei.posB + 1)
        loci_similar = False
        loci_same = False
        intra = ei.chrA == ei.chrB == ej.chrA == ej.chrB
        if intra:
            loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posA, ej.posB, ei.preciseA, ei.preciseB, ej.preciseA, ej.preciseB, min(ei.svlen, ej.svlen))
        elif ei.chrA == ej.chrA and ei.chrB == ej.chrB:  # Try chrA matches chrA
            loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posA, ej.posB, ei.preciseA, ei.preciseB, ej.preciseA, ej.preciseB, min(ei.svlen, ej.svlen), intra=False)
        elif ei.chrA == ej.chrB and ei.chrB == ej.chrA:
            loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posB, ej.posA, ei.preciseA, ei.preciseB, ej.preciseB, ej.preciseA, min(ei.svlen, ej.svlen), intra=False)
        if not loci_similar:
            continue
        # Force merging of translocations that have similar loci
        if not intra:
            G.add_edge(i_id, j_id, loci_same=loci_same)  #, w=0)
            continue

        one_is_imprecise = False
        both_imprecise = False
        if (not ei.preciseA or not ei.preciseB) and ei.svlen_precise:
            one_is_imprecise = True  # also not remapped
        if (not ej.preciseA or not ej.preciseB) and ej.svlen_precise:
            if one_is_imprecise:
                both_imprecise = True
            else:
                one_is_imprecise = True
        if paired_end and ei.su == ej.su == 1 and not ej.preciseA and not ej.preciseB:
            continue

        ci = ei.contig
        ci2 = ei.contig2
        ci_alt = ei.variant_seq
        if isinstance(ci_alt, str):
            if len(ci_alt) == 1 or (len(ci_alt) > 0 and ci_alt[0] in "<."):
                ci_alt = ""
        cj = ej.contig
        cj2 = ej.contig2
        cj_alt = ej.variant_seq
        if isinstance(cj_alt, str):
            if len(cj_alt) == 1 or (len(cj_alt) > 0 and cj_alt[0] in "<."):
                cj_alt = ""

        any_contigs_to_check = any((ci, ci2, ci_alt)) and any((cj, cj2, cj_alt))

        # dont merge deletion events with low/poor evidence
        if same_sample and not loci_same and ei.svtype == "DEL" and ei.su < 3 and ej.su < 3 and not any_contigs_to_check and ei.spanning == 0 and ej.spanning == 0 and ei.sc == 0 and ej.sc == 0:
            continue
        recpi_overlap = is_reciprocal_overlapping(ei.posA, ei.posB, ej.posA, ej.posB)
        overlap = max(0, min(ei.posA, ej.posA) - max(ei.posB, ej.posB))
        if paired_end:
            if ei.spanning > 0 and ej.spanning > 0 and overlap == 0:
                disjoint_nodes.add(i_id)
                disjoint_nodes.add(j_id)
                continue
        m = False
        ml = max(int(ei.svlen), int(ej.svlen))
        if ml == 0 and ei.svtype != 'TRA':
            continue
        if ei.svtype == 'TRA':
            l_ratio = 1
        else:
            l_ratio = min(int(ei.svlen), int(ej.svlen)) / ml
        if ei.svtype == "INS":
            if aggressive_ins_merge:
                m = True
            else:
                # if both have been remapped, make sure size is similar
                if ei.spanning > 0 or ej.spanning > 0:
                    if paired_end:
                        if isinstance(ei.variant_seq, str) and isinstance(ej.variant_seq, str):
                            if l_ratio > 0.9:
                                m = True
                        else:
                            m = True
                    elif ml > 0 and l_ratio > 0.8:
                        m = True
                elif ei.remap_score > 0 and ej.remap_score > 0 and ml > 0 and l_ratio > 0.8:
                    m = True
                elif ei.remap_score == 0 or ej.remap_score == 0:
                    m = True
        elif ei.svtype != "INS":
            spd = span_position_distance(ei, ej)
            if paired_end or aggressive_ins_merge:  # when aggresive_ins_merge is True, then this is dysgu merge, not dysgu call
                if ei.spanning > 0 and ej.spanning > 0:
                    if recpi_overlap or spd > 0.3:
                        m = True
                elif (recpi_overlap or spd > 0.3 or (loci_similar and any_contigs_to_check)) and not both_in_include:
                    m = True
            elif (recpi_overlap or spd > 0.3 or (loci_similar and any_contigs_to_check)) and not both_in_include:
                m = True
        if not m:
            continue
        # Loci are similar, check contig match or reciprocal overlap
        if not any_contigs_to_check:
            if ml > 0:
                if l_ratio > 0.5 or (one_is_imprecise and l_ratio > 0.3):
                    G.add_edge(i_id, j_id, loci_same=loci_same)
        else:
            v = None
            if ci_alt and cj_alt:
                v = (ci_alt, cj_alt)
            elif ci and cj:
                v = (ci, cj)
            elif ci and cj2:
                v = (ci, cj2)
            elif ci2 and cj2:
                v = (ci2, cj2)
            elif ci2 and cj:
                v = (ci2, cj)
            else:
                continue
            # echo(ei.remap_score == 0 or ej.remap_score == 0, ei.svlen, ej.svlen, v, assembler.check_contig_match(v[0], v[1], return_int=True))
            # if not remapped and nearby insertions with opposing soft-clips --> merge
            # also long-reads will normally have remap_score == 0
            if ei.remap_score == 0 or ej.remap_score == 0:
                if (v[0][0].islower() and v[1][-1].islower()) or (v[0][-1].islower() and v[1][0].islower()):
                    G.add_edge(i_id, j_id, loci_same=loci_same)
                    continue
            if assembler.check_contig_match(v[0], v[1], return_int=True):
                G.add_edge(i_id, j_id, loci_same=True)
            # see if alt sequence can be found in other contig
            else:
                if ci_alt and cj:
                    if assembler.check_contig_match(ci_alt, cj, return_int=True):
                        G.add_edge(i_id, j_id, loci_same=True)
                        continue
                if ci_alt and cj2:
                    if assembler.check_contig_match(ci_alt, cj2, return_int=True):
                        G.add_edge(i_id, j_id, loci_same=True)
                        continue
                if cj_alt and ci:
                    if assembler.check_contig_match(cj_alt, ci, return_int=True):
                        G.add_edge(i_id, j_id, loci_same=True)
                        continue
                if cj_alt and ci2:
                    if assembler.check_contig_match(cj_alt, ci2, return_int=True):
                        G.add_edge(i_id, j_id, loci_same=True)
                        continue
    return G, disjoint_nodes


def cut_components(G, disjoint_nodes):
    # e = G.edges(data=True)
    # G2 = nx.Graph([i for i in e if i[2]["loci_same"] == True])
    # for u in G.nodes():
    #     if u not in G2:
    #         e0 = next(G.edges(u).__iter__())  # Use first edge out of u to connect
    #         G2.add_edge(*e0)
    # components = nx.algorithms.components.connected_components(G2)
    components = nx.algorithms.components.connected_components(G)
    if len(disjoint_nodes) > 0:
        # try split this component into disjoint sets. This method works for small cluster sizes (most of the time)
        # but can fail when there are many disjoint nodes. Label propagation might be needed for these
        components2 = []
        for c in components:
            n_disjoin = set([])
            for node in c:
                if node in disjoint_nodes:
                    n_disjoin.add(node)
            if len(n_disjoin) <= 1:
                components2.append(c)
                continue
            out_e = defaultdict(list)
            for node in n_disjoin:
                for neigh in G.neighbors(node):
                    out_e[neigh].append(node)
            G3 = nx.Graph()
            for k, v in out_e.items():
                G3.add_edge(k, random.choice(v))
            components2 += list(nx.algorithms.components.connected_components(G3))
        return components2
    return components


cpdef srt_func(c):
    if c.type != "pe" and c.type != "":
        return 100 + c.su
    return c.su + (3 * c.spanning)


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
    G, disjoint_nodes = enumerate_events(G, potential, max_dist, try_rev, tree, paired_end, rel_diffs, diffs, same_sample,
                        aggressive_ins_merge=aggressive_ins_merge,
                        debug=debug)
    found = []
    for item in potential:  # Add singletons, non-merged
        if not G.has_node(item.event_id):
            found.append(item)
    # Try and merge SVs with identical breaks, then merge ones with less accurate breaks - this helps prevent
    # over merging SVs that are close together
    components = cut_components(G, disjoint_nodes)
    node_to_event = {i.event_id: i for i in potential}
    cdef int k
    for grp in components:
        c = [node_to_event[n] for n in grp]
        best = sorted(c, key=srt_func, reverse=True)
        w0 = best[0]
        if not pick_best:
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
            sqc = w0.sqc
            best_var_seq = w0.variant_seq
            for k in range(1, len(best)):
                item = best[k]
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
                        elif min_size > w0.svlen < item.svlen:
                            w0.svlen = item.svlen
                elif item.svtype == "INS" and svt == "INS":
                    if not spanned:
                        if item.spanning:
                            w0.svlen = item.svlen
                            w0.variant_seq = item.variant_seq
                        elif item.svlen * 0.6 < w0.svlen < item.svlen or min_size > w0.svlen < item.svlen:
                            if best_var_seq == -1:
                                w0.svlen = item.svlen
                        if best_var_seq is None and isinstance(item.variant_seq, str):
                            w0.variant_seq = item.variant_seq
                            w0.svlen = item.svlen
                            best_var_seq = item.variant_seq
                        if w0.posA == w0.posB < item.posB:
                            w0.posB = item.posB
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
                    w0.MAPQsupp = norm_vals(w0.MAPQsupp, weight, item.MAPQsupp, wt, denom)
                    w0.MAPQpri = norm_vals(w0.MAPQpri, weight, item.MAPQpri, wt, denom)
                    w0.NMpri = norm_vals(w0.NMpri, weight, item.NMpri, wt, denom)
                    w0.NMsupp = norm_vals(w0.NMsupp, weight, item.NMsupp, wt, denom)
                if add_contig_a and item.contig and len(item.contig) > len(new_a):
                    new_a = item.contig
                if add_contig_b and item.contig2 and len(item.contig2) > len(new_b):
                    new_b = item.contig2
            if add_contig_a and new_a:
                w0.contig = new_a
            if add_contig_b and new_b:
                w0.contig2 = new_b
        if add_partners:
            w0.partners = [i.event_id for i in best[1:]]
        found.append(w0)
    return found


def sample_level_density(potential, regions, max_dist=50):
    tmp_list = defaultdict(list)
    cdef EventResult_t ei
    for idx in range(len(potential)):
        ei = potential[idx]
        # Only find density for non-region calls, otherwise too dense to be meaningful
        if not intersecter(regions, ei.chrA, ei.posA, ei.posA + 1):
            tmp_list[ei.chrA].append((ei.posA - max_dist, ei.posA + max_dist, idx))
        if not intersecter(regions, ei.chrB, ei.posB, ei.posB + 1):
            tmp_list[ei.chrB].append((ei.posB - max_dist, ei.posB + max_dist, idx))
    nc2 = {k: io_funcs.iitree(v, add_value=True) for k, v in tmp_list.items()}
    cdef int vv
    for idx in range(len(potential)):
        ei = potential[idx]
        neighbors = 0.
        count = 0.
        if ei.chrA == ei.chrB and abs(ei.posB - ei.posA) < 2:
            expected = 2
        else:
            expected = 1
        if not intersecter(regions, ei.chrA, ei.posA, ei.posA + 1):
            vv = nc2[ei.chrA].countOverlappingIntervals(ei.posA, ei.posA + 1)
            neighbors += vv - expected
            count += 1
        if not intersecter(regions, ei.chrB, ei.posB, ei.posB + 1):
            vv = nc2[ei.chrB].countOverlappingIntervals(ei.posB, ei.posB + 1)
            neighbors += vv - expected
            count += 1
        neighbors_10kb = 0.
        count_10kb = 0
        large_itv = merge_intervals(((ei.chrA, ei.posA, ei.posA + 1), (ei.chrB, ei.posB, ei.posB + 1)), pad=10000)
        for c, s, e in large_itv:
            if not intersecter(regions, c, s, e):
                vv = nc2[c].countOverlappingIntervals(s, e)
                neighbors_10kb += vv - len(large_itv)
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
                for j in range(i, finish, t):
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
                  merge_dist, regions_only, assemble_contigs, rel_diffs, diffs, min_size, max_single_size,
                  sites_index, paired_end, length_extend, divergence):
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
                                                      divergence):
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
        max_single_size, sites_index, paired_end, length_extend, divergence = args
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
                                                   paired_end=paired_end, length_extend=length_extend, divergence=divergence)
        if potential_events:
            for res in potential_events:
                pickle.dump(res, completed_file)
    completed_file.close()


# def postcall_job(preliminaries, aux_data):
#     mode, ref_path, no_gt, insert_stdev, paired_end, drop_gaps = aux_data
#     ref_genome = pysam.FastaFile(ref_path)
#     preliminaries = re_map.drop_svs_near_reference_gaps(preliminaries, paired_end, ref_genome, drop_gaps)
#     preliminaries = post_call.ref_repetitiveness(preliminaries, ref_genome)
#     preliminaries = post_call.strand_binom_t(preliminaries)
#     preliminaries = assembler.contig_info(preliminaries)  # GC info, repetitiveness
#     preliminaries = find_repeat_expansions(preliminaries, insert_stdev)
#     preliminaries = post_call.compressability(preliminaries)
#     preliminaries = post_call.get_gt_metric2(preliminaries, mode, no_gt)
#     preliminaries = post_call.get_ref_base(preliminaries, ref_genome)
#     return preliminaries


# https://rvprasad.medium.com/data-and-chunk-sizes-matter-when-using-multiprocessing-pool-map-in-python-5023c96875ef
# aux_data = None
#
# def initializer(init_data):
#     global aux_data
#     aux_data = init_data


# def with_initializer_worker_wrapper(varying_data):
#     return postcall_job(varying_data, aux_data)


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
        logging.info(f"Max clustering dist {max_clust_dist}")
        args["divergence"] = 1
    else:
        if args["divergence"] == "auto":
            divergence, _ = genome_scanner.get_read_properties(args["max_tlen"], insert_median, insert_stdev, read_len, ibam, find_divergence=True)
            args["divergence"] = divergence
        else:
            args["divergence"] = float(args["divergence"])
            logging.info(f"Sequence divergence upper bound {args['divergence']}")
        if args["mode"] == "pacbio":
            max_dist, max_clust_dist = 35, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 700
        elif args["mode"] == "nanopore":
            max_dist, max_clust_dist = 100, 500000
            if args["merge_dist"] is None:
                args["merge_dist"] = 700

    # set upper bound on single-partition size
    max_single_size = min(max(args["max_cov"] * 50, 10000), 100000)  # limited between 5000 - 50,000 reads
    event_id = 0
    block_edge_events = []

    clip_length = args["clip_length"]
    merge_dist = args["merge_dist"]
    min_size = args["min_size"]
    length_extend = args["length_extend"]
    divergence = args["divergence"]
    read_buffer = genome_scanner.read_buffer
    sites_info = sites_utils.vcf_reader(args["sites"], infile, args["parse_probs"], sample_name, args["ignore_sample_sites"] == "True", args["sites_prob"], args["sites_pass_only"] == "True")

    cdef Py_SimpleGraph G
    G, node_to_name, bad_clip_counter, sites_adder, n_aligned_bases = graph.construct_graph(genome_scanner,
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
                                            find_n_aligned_bases=find_n_aligned_bases)
    sites_index = None
    if sites_adder:
        sites_index = sites_adder.sites_index
    logging.info("Graph constructed")

    auto_support = False
    if args["min_support"] != "auto":
        args["min_support"] = int(args["min_support"])
        min_support = args["min_support"]
        logging.info(f"Minimum support {args['min_support']}")
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
        logging.info(f"Inferred minimum support {min_support}")

    if args["pl"] == "pe":  # reads with internal SVs can be detected at lower support
        lower_bound_support = min_support - 1 if min_support - 1 > 1 else 1
    else:
        lower_bound_support = min_support

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
                rel_diffs, diffs, min_size, max_single_size, sites_index, paired_end, length_extend, divergence )
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
                                                                   divergence=divergence)
                        if potential_events:
                            if not low_mem:
                                block_edge_events += potential_events
                            else:
                                for res in potential_events:
                                    pickle.dump(res, completed_file)
                    else:
                        j_submitted, w_idx = heapq.heappop(minhq)
                        heapq.heappush(minhq, (j_submitted + len(res["n2n"]), w_idx))
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
                                                                   length_extend=length_extend, divergence=divergence)
                        if potential_events:
                            if not low_mem:
                                block_edge_events += potential_events
                            else:
                                for res in potential_events:
                                    pickle.dump(res, completed_file)
                    else:
                        j_submitted, w_idx = heapq.heappop(minhq)
                        heapq.heappush(minhq, (j_submitted + len(res["n2n"]), w_idx))
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
    logging.info("Number of components {}. N candidates {}".format(components_seen, len(block_edge_events)))
    keeps = len([i for i in block_edge_events if i.site_info])
    if keeps:
        logging.info("Number of matching SVs from --sites {}".format(keeps))
    preliminaries = []
    if args["remap"] == "True" and args["contigs"] == "True":
        block_edge_events = re_map.remap_soft_clips(block_edge_events, ref_genome,
                                                    keep_unmapped=True if args["pl"] == "pe" else False,
                                                    min_support=min_support)
        logging.info("Re-alignment of soft-clips done. N candidates {}".format(len(block_edge_events)))
    # Merge across calls
    if args["merge_within"] == "True":
        merged = merge_events(block_edge_events, args["merge_dist"], regions, bool(paired_end), try_rev=False, pick_best=False,
                                         debug=True, min_size=args["min_size"])
    else:
        merged = block_edge_events
    logging.info("Number of candidate SVs merged: {}".format(len(block_edge_events) - len(merged)))
    logging.info("Number of candidate SVs after merge: {}".format(len(merged)))
    before = len(merged)

    # if not args["keep_small"]:
    #     merged = [event for event in merged if (event.svlen >= args["min_size"] or event.chrA != event.chrB) and (event.su >= args["min_support"] or event.site_info)]
    #     logging.info("Number of candidate SVs dropped with sv-len < min-size or support < min support: {}".format(before - len(merged)))
    # else:
    #     merged = [event for event in merged if (event.su >= args["min_support"] or event.site_info)]
    #     logging.info("Number of candidate SVs dropped with support < min support: {}".format(before - len(merged)))
    #
    # cdef EventResult_t d
    # for d in merged:
    #     d.type = args["pl"]
    # merged = re_map.drop_svs_near_reference_gaps(merged, paired_end, ref_genome, args["drop_gaps"] == "True")
    # coverage_analyser = post_call.CoverageAnalyser(tdir)
    #
    # preliminaries = coverage_analyser.process_events(merged)
    # preliminaries = coverage.get_raw_coverage_information(merged, regions, coverage_analyser, infile, args["max_cov"])
    # preliminaries = sample_level_density(preliminaries, regions)
    # preliminaries = post_call.get_badclip_metric(preliminaries, bad_clip_counter, infile, regions)
    # preliminaries = post_call.ref_repetitiveness(preliminaries, ref_genome)
    # preliminaries = post_call.strand_binom_t(preliminaries)
    # preliminaries = assembler.contig_info(preliminaries)  # GC info, repetitiveness
    # preliminaries = find_repeat_expansions(preliminaries, insert_stdev)
    # preliminaries = post_call.compressability(preliminaries)
    # preliminaries = post_call.get_gt_metric2(preliminaries, args["mode"], args["no_gt"])
    # preliminaries = coverage_analyser.normalize_coverage_values(preliminaries)
    # preliminaries = post_call.get_ref_base(preliminaries, ref_genome)



    if auto_support:
        if not args["keep_small"]:
            merged = [event for event in merged if (event.svlen >= args["min_size"] or event.chrA != event.chrB) or event.site_info]
    else:
        if not args["keep_small"]:
            merged = [event for event in merged if (event.svlen >= args["min_size"] or event.chrA != event.chrB) and (event.su >= args["min_support"] or event.site_info)]
        else:
            merged = [event for event in merged if (event.su >= args["min_support"] or event.site_info)]

    cdef EventResult_t d
    for d in merged:
        d.type = args["pl"]

    coverage_analyser = post_call.CoverageAnalyser(tdir)
    preliminaries = coverage_analyser.process_events(merged)
    preliminaries = coverage.get_raw_coverage_information(merged, regions, coverage_analyser, infile, args["max_cov"])

    if auto_support:
        preliminaries = post_call.filter_auto_min_support(preliminaries)
        logging.info("Number of candidate SVs dropped with support < min support: {}".format(before - len(merged)))

    preliminaries = post_call.get_badclip_metric(preliminaries, bad_clip_counter, infile, regions)

    # if False:
    #     aux_data = args["mode"], args["reference"], args["no_gt"], insert_stdev, paired_end, args["drop_gaps"] == "True"
    #     pool = multiprocessing.Pool(procs, initializer, (aux_data,))
    #     n = int(len(preliminaries) / procs)
    #     preliminaries = [preliminaries[i:i + n] for i in range(0, len(preliminaries), n)]
    #     preliminaries = pool.map(with_initializer_worker_wrapper, preliminaries)
    #     preliminaries = [item for m in preliminaries for item in m]
    #
    # else:

    preliminaries = re_map.drop_svs_near_reference_gaps(preliminaries, paired_end, ref_genome, args["drop_gaps"] == "True")
    preliminaries = post_call.ref_repetitiveness(preliminaries, ref_genome)
    preliminaries = post_call.strand_binom_t(preliminaries)
    preliminaries = assembler.contig_info(preliminaries)  # GC info, repetitiveness
    preliminaries = find_repeat_expansions(preliminaries, insert_stdev)
    preliminaries = post_call.compressability(preliminaries)
    preliminaries = post_call.get_gt_metric2(preliminaries, args["mode"], args["no_gt"])
    preliminaries = post_call.get_ref_base(preliminaries, ref_genome, args["symbolic_sv_size"])

    preliminaries = sample_level_density(preliminaries, regions)
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

    #####################
    #  Run dysgu here   #
    #####################
    events, site_adder = pipe1(args, infile, kind, regions, ibam, ref_genome, sample_name)
    if not events:
        logging.critical("No events found")
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
            contig_header_lines = ""
            for item in infile.header["SQ"]:
                contig_header_lines += f"\n##contig=<ID={item['SN']},length={item['LN']}>"
            args["add_kind"] = "True"
            args["sample_name"] = sample_name
            io_funcs.to_vcf(df, args, {sample_name}, outfile, show_names=False, contig_names=contig_header_lines,
                            sort_output=False)
    logging.info("dysgu call {} complete, n={}, time={} h:m:s".format(
               args["sv_aligns"],
               len(df),
               str(datetime.timedelta(seconds=int(time.time() - t0)))))
