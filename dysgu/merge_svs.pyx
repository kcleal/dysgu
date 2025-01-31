# cython: language_level=3

from __future__ import absolute_import

from multiprocessing import Pool

import numpy as np
import random
from collections import defaultdict, deque
import networkx as nx
from dysgu import consensus
from dysgu.map_set_utils cimport is_reciprocal_overlapping, EventResult
from dysgu.map_set_utils import echo
from cython.operator import dereference
from functools import cmp_to_key
from superintervals import IntervalSet


ctypedef EventResult EventResult_t

np.random.seed(0)
random.seed(0)

def get_chrom_key(ei):
    if ei.chrA != ei.chrB:
        return min(ei.chrA, ei.chrB), max(ei.chrA, ei.chrB)
    return ei.chrA


def compare_subset(potential, max_dist, max_comparisons, same_sample):
    tmp_list = defaultdict(list)
    cdef int idx, jdx, dist1, ci_a, start, stop
    cdef int half_d = <int>(max_dist * 0.5)
    for idx in range(len(potential)):
        ei = potential[idx]
        chrom_key = get_chrom_key(ei)
        tmp_list[chrom_key].append((ei.posA - ei.cipos95A - max_dist, ei.posA + ei.cipos95A + max_dist, idx))
        if ei.chrA == ei.chrB and ei.svlen > half_d:
            tmp_list[chrom_key].append(
                (ei.posB - ei.cipos95B - max_dist, ei.posB + ei.cipos95B + max_dist, idx))

    si_sets = {}
    for k, v in tmp_list.items():
        iset = IntervalSet(with_data=True)
        for start, stop, idx in v:
            iset.add_int_value(start, stop, idx)
        iset.index()
        si_sets[k] = iset

    for idx in range(len(potential)):
        ei = potential[idx]
        # ols = nc2[get_chrom_key(ei)].allOverlappingIntervals(ei.posA, ei.posA + 1)
        ols = si_sets[get_chrom_key(ei)].find_overlaps(ei.posA, ei.posA + 1)
        # echo(ols, ols2)
        ols = [i for i in set(ols) if i != idx]
        if len(ols) > max_comparisons:
            random.shuffle(ols)
            ols = ols[:max_comparisons]
        for jdx in ols:
            ej = potential[jdx]
            yield ei, ej, idx, jdx


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
    cdef int temp, dist_a, dist_b, precise_thresh, imprecise_thresh, same_thresh
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


cpdef similar_locations(intra, ei, ej):
    if intra:
        loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posA, ej.posB, ei.preciseA, ei.preciseB,
                                                  ej.preciseA, ej.preciseB, min(ei.svlen, ej.svlen))
    elif ei.chrA == ej.chrA and ei.chrB == ej.chrB:  # Try chrA matches chrA
        loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posA, ej.posB, ei.preciseA, ei.preciseB,
                                                  ej.preciseA, ej.preciseB, min(ei.svlen, ej.svlen), intra=False)
    elif ei.chrA == ej.chrB and ei.chrB == ej.chrA:
        loci_similar, loci_same = break_distances(ei.posA, ei.posB, ej.posB, ej.posA, ei.preciseA, ei.preciseB,
                                                  ej.preciseB, ej.preciseA, min(ei.svlen, ej.svlen), intra=False)
    if not loci_similar:
        return False
    return loci_similar, loci_same


def consistent_alignment_and_cigars(ei, ej, l_ratio):
    # if contig location is set, check for overlapping alignment
    if ei.contig_ref_end > 0 and ej.contig_ref_end > 0:
        ref_overlap = min(ei.contig_ref_end, ej.contig_ref_end) - max(ei.contig_ref_start, ej.contig_ref_start)
        if ref_overlap < 0:
            return False
    else:
        return True

    # skip merging for very repetitive seqs
    if ei.rep > 0.3 and ej.rep > 0.3 and l_ratio < 0.5:
        return False

    # skip merging if contigs support multiple events or very different events
    if not ei.contig_cigar or not ej.contig_cigar:
        return True
    if ej.svtype == "INS":
        diff_i = sum(l for op, l in ei.contig_cigar if op == 1 and l > 10)
        diff_j = sum(l for op, l in ej.contig_cigar if op == 1 and l > 10)
    if ej.svtype == "DEL":
        diff_i = sum(l for op, l in ei.contig_cigar if op == 2 and l > 10)
        diff_j = sum(l for op, l in ej.contig_cigar if op == 2 and l > 10)
    tot_non_ref = diff_i + diff_j
    tot_sv_len = ei.svlen + ej.svlen
    ratio = tot_non_ref / tot_sv_len
    if ratio > 1.8: # or ratio < 0.25:
        return False
    return True


cdef float jaccard_similarity(set1, set2):
    if not set1 or not set2:
        return 0
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0


cpdef get_consensus_seqs(ei, ej):
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
    return any_contigs_to_check, ci, ci2, ci_alt, cj, cj2, cj_alt


def contig_pairs_iter(ci, ci2, ci_alt, cj, cj2, cj_alt):
    if ci and cj:
        yield ci, cj
    if ci2 and cj2:
        yield ci2, cj2
    if ci2 and cj:
        yield ci2, cj
    if ci and cj2:
        yield ci, cj2
    if ci_alt and cj:
        yield ci_alt, cj
    if ci_alt and cj2:
        yield ci_alt, cj2
    if ci and cj_alt:
        yield ci, cj_alt
    if ci2 and cj_alt:
        yield ci2, cj_alt
    if ci_alt and cj_alt:
        yield ci_alt, cj_alt


cdef int matches_with_max_gap(char *c_cigar, int max_gap, float gaps_tol ):
    # get matching bases, and filter large alignment gaps
    cdef:
        int matches = 0
        int gaps = 0
        int num = 0  # To accumulate the number as we parse through the string
        char op

    while dereference(c_cigar):
        if b'0' <= dereference(c_cigar) <= b'9':
            # Convert digit char to int and accumulate
            num = num * 10 + (dereference(c_cigar) - 48 )
        else:
            op = dereference(c_cigar)
            if op == b'D' or op == b'I':
                if num > max_gap:
                    # echo('max gap', num)
                    return 0
                gaps += num
            if op == b'M':
                matches += num
            num = 0
        c_cigar += 1
    if matches and <float>gaps / <float>matches > gaps_tol:
        # echo('gaps vs matches', gaps / matches)
        return 0
    return matches


cdef bint bad_insertion_coverage(char *c_cigar, seq, bint is_query, float ins_gap_tol):
    # Count the matching bases over the lowercase 'insertion' seq only
    cdef:
        int seq_index = 0
        int matches = 0
        int gaps = 0
        int num = 0
        int i
        char op

    while dereference(c_cigar):
        if b'0' <= dereference(c_cigar) <= b'9':
            num = num * 10 + (dereference(c_cigar) - 48)
        else:
            op = dereference(c_cigar)
            if op == b'D' or op == b'I':
                if seq[seq_index].islower():
                    gaps += num
                if op == b'D' and not is_query:
                    seq_index += num
            if op == b'M':
                for i in range(seq_index, seq_index + num - 1):
                    if seq[i].islower():
                        matches += 1
                seq_index += num
            num = 0
        c_cigar += 1
    if not matches or <float> gaps / <float> matches > ins_gap_tol:
        # echo('ins cov', gaps, matches)
        return <bint>True
    return <bint>False


cdef bint bad_alignment(alignment, ei, ej, v, paired_end,
                        float soft_clip_tol, float ins_gap_tol, float gaps_tol,
                        int max_gap):
    # Check if alignment is consistent. ej is optional
    if not alignment:
        return True

    qs = alignment.query_begin
    qe = alignment.query_end
    qlen = len(v[0])
    tlen = len(v[1])
    ts = alignment.target_begin
    te = alignment.target_end_optimal

    q_right_clip = qlen - qe
    t_right_clip = tlen - te

    if not paired_end:
        soft_clip_tolerance = max(10, min(qlen, tlen) * soft_clip_tol)
        matches_threshold = min(qlen, tlen) * 0.6
        ins_cov = 0.8
    else:
        soft_clip_tolerance = 10
        matches_threshold = min(qlen, tlen) * 0.5
        ins_cov = 0.9

    cdef bytes t = alignment.cigar.encode('utf8')
    cdef char *c_cigar = t

    # Skip alignments with larger gaps or not many matching bases
    cdef int matches = matches_with_max_gap(c_cigar, max_gap, gaps_tol)
    if matches == 0 or matches < matches_threshold:
        # echo('no match threshold', matches)
        return True

    # Keep if aligned inserted bases ~ svlen
    if (ei.svtype == "INS" and ei.spanning > 0) or (ej and ej.svtype == "INS" and ej.spanning > 0):
        a_into_b = True
        b_into_a = True
        if ei.svtype == "INS" and ei.spanning > 0:
            ei_ins_count = len([i for i in v[0][qs:qe] if i.islower()])
            if ei_ins_count / ei.svlen < ins_cov:
                a_into_b = False
            elif bad_insertion_coverage(c_cigar, v[0][qs:qe], <bint>True, ins_gap_tol):
                a_into_b = False
        if ej:
            if ej.svtype == "INS" and ej.spanning > 0:
                ej_ins_count = len([i for i in v[1][ts:te] if i.islower()])
                if bad_insertion_coverage(c_cigar, v[1][ts:te], <bint>False, ins_gap_tol):
                    b_into_a = False
                elif ej_ins_count / ej.svlen < ins_cov:
                    b_into_a = False
            if not a_into_b and not b_into_a:
                # echo('not a into b')
                return True

    # Keep, if one sequence is more or less completely aligned to the other
    if qs < soft_clip_tolerance and q_right_clip < soft_clip_tolerance:
        return False
    if ts < soft_clip_tolerance and t_right_clip < soft_clip_tolerance:
        return False

    if qs > soft_clip_tolerance:
        if q_right_clip > soft_clip_tolerance:
            # echo('tol 1')
            return True
        if ts > soft_clip_tolerance:
            # echo('tol 1b', ts, soft_clip_tolerance)
            return True

    if ts > soft_clip_tolerance:
        if t_right_clip > soft_clip_tolerance:
            # echo('tol 3')
            return True
        if qs > soft_clip_tolerance:
            # echo('tol 4')
            return True

    if q_right_clip > soft_clip_tolerance:
        if qs > soft_clip_tolerance:
            # echo('tol 5')
            return True
        if t_right_clip > soft_clip_tolerance:
            # echo('tol 6')
            return True

    if t_right_clip > soft_clip_tolerance:
        if ts > soft_clip_tolerance:
            # echo('tol 7')
            return True
        if q_right_clip > soft_clip_tolerance:
            # echo('tol 8')
            return True

    return False  # ok


def process_contig_aignments(ci, ci2, ci_alt, cj, cj2, cj_alt, ei, ej, paired_end, idx, jdx, same_sample):
    cdef float soft_clip_tol, ins_gaps_tol, gaps_tol
    cdef int max_gap
    if ei.type == "nanopore" or ej.type == "nanopore":
        soft_clip_tol = 0.05
        ins_gaps_tol = 0.08
        gaps_tol = 0.07
        max_gap = 14
    else:
        soft_clip_tol = 0.04
        ins_gaps_tol = 0.05
        gaps_tol = 0.05
        max_gap = 12

    for v in contig_pairs_iter(ci, ci2, ci_alt, cj, cj2, cj_alt):
        if not v[0] or not v[1]:
            continue
        if same_sample:
            res = consensus.check_contig_match(v[0], v[1], return_int=False)
            # echo(res)
            if bad_alignment(res, ei, ej, v, paired_end, soft_clip_tol, ins_gaps_tol, gaps_tol, max_gap):
                # echo('bad alignment')
                break
            # echo("---->MERGED3", ei.svlen, ej.svlen, ei.svtype, (idx, jdx))
            return idx, jdx

        elif consensus.check_contig_match(v[0], v[1], return_int=True):
            return idx, jdx


# def consistent_break_pattern(ei, ej, ci, ci2, cj, cj2):  # for paired-end reads only
#     return True
#     if ei.svtype == ej.svtype:
#         return True
#     if ei.svtype in 'INSTRABND' or ei.svtype in 'INSTRABND':
#         if ei.svtype in 'INSTRABND':
#             if ej.svtype == 'DUP':  # Make sure merging SV has clip on left-hand-side
#                 if ej.posA <= ei.posA and ci and ci[0].islower():
#                     return True
#             if ej.svtype == 'DEL':
#                 if ej.posA >= ei.posA and ci and ci[-1].islower():
#                     return True
#
#         elif ej.svtype in 'INSTRABND':
#             if ei.svtype == 'DUP':
#                 if ei.posA <= ej.posA and cj and cj[0].islower():
#                     return True
#             if ei.svtype == 'DEL':
#                 if ei.posA >= ej.posA and cj and cj[-1].islower():
#                     return True
#     return False  # No merging between other SV types


def enumerate_events(G, potential, max_dist, try_rev, tree, paired_end=False, rel_diffs=False, diffs=15,
                     same_sample=True, aggressive_ins_merge=False, debug=False, max_comparisons=20, procs=1):
    event_iter = compare_subset(potential, max_dist, max_comparisons, same_sample)

    seen, disjoint_edges = set(), set()
    out_edges = defaultdict(int)

    avg_su_thresh = max(np.median([i.su for i in potential]) * 0.4, 10)

    job = []
    pool = None
    if procs > 1:
        pool = Pool(procs)
    for ei, ej, idx, jdx in event_iter:

        if len(job) > 1000:
            if procs == 1:
                for item in job:
                    edge = process_contig_aignments(*item)
                    if not edge:
                        continue
                    else:
                        out_edges[edge[0]] += 1
                        out_edges[edge[1]] += 1
                        G.add_edge(edge[0], edge[1], loci_same=True)
            else:
                for edge in pool.starmap(process_contig_aignments, job):
                    if not edge:
                        continue
                    else:
                        out_edges[edge[0]] += 1
                        out_edges[edge[1]] += 1
                        G.add_edge(edge[0], edge[1], loci_same=True)
            job = []

        i_id, j_id = ei.event_id, ej.event_id

        id_key = (min(idx, jdx), max(idx, jdx))
        idx, jdx = id_key

        if id_key in seen:
            continue

        if not same_sample and ((ei.svtype in {"INS", "DUP"} and ej.svtype not in {"INS", "DUP"}) or ei.svtype != ej.svtype):
            seen.add(id_key)
            continue
        #todo
        if same_sample and \
            ei.su > avg_su_thresh and ej.su > avg_su_thresh and \
            ei.spanning > 0 and ej.spanning > 0 and \
            ei.svlen != ej.svlen:  # <-- todo THIS
            seen.add(id_key)
            disjoint_edges.add(id_key)
            # echo('disjoin added here', ei.su, ej.su, avg_su_thresh)
            continue

        if not same_sample and ei.sample == ej.sample:
            seen.add(id_key)
            disjoint_edges.add(id_key)
            continue

        if (same_sample and i_id == j_id) or \
            (not same_sample and ei.sample == ej.sample) or \
            out_edges[idx] > 10 or out_edges[jdx] > 10:
                # or
                # (ej.svtype == "DEL" and ei.svtype != "DEL") or
                # (ei.svtype == "DEL" and ej.svtype != "DEL")):
                # (ei.svtype in {"INS", "DUP"} and ej.svtype not in {"INS", "DUP"}) or
                # ei.svtype != ej.svtype):
                seen.add(id_key)
                continue

        seen.add(id_key)


        # if ei.posA == 136934:
        # if id_key not in seen:
        # echo("merge candidate", (id_key), ei.svlen, ej.svlen, "positions", ei.posA, ej.posA, ei.svtype, ej.svtype, 'SU=', ei.su, ej.su, ei.haplotypes, ej.haplotypes)


        intra = ei.chrA == ej.chrA and ei.chrB == ej.chrB

        # loci_similar = True
        # if paired_end:
        #     loci_similar = similar_locations(intra, ei, ej)
        #     if not loci_similar:
        #         continue

        if not intra:
            out_edges[idx] += 1
            out_edges[jdx] += 1
            G.add_edge(idx, jdx)
            continue

        if ei.type != ej.type and (ei.type == 'pe' or ej.type == 'pe'):
            out_edges[idx] += 1
            out_edges[jdx] += 1
            G.add_edge(idx, jdx)
            continue

        if ei.spanning > 0 and ej.spanning > 0 and jaccard_similarity(ei.qnames, ej.qnames) > 0.1: #0.1:
            # echo('qnames!', ei.qnames, ej.qnames, jaccard_similarity(ei.qnames, ej.qnames))
            disjoint_edges.add(id_key)
            continue

        if ei.posA == ej.posA and ei.svtype == ej.svtype and ei.svlen == ej.svlen:
            out_edges[idx] += 1
            out_edges[jdx] += 1
            G.add_edge(idx, jdx, loci_same=True)

        any_contigs_to_check, ci, ci2, ci_alt, cj, cj2, cj_alt = get_consensus_seqs(ei, ej)

        if paired_end:

            # Make sure patterns are consistent between different SV calls
            # if not consistent_break_pattern(ei, ej, ci, ci2, cj, cj2):
            #     continue
            # if paired_end and ei.svtype != ej.svtype:
            #     echo("merge candidate", (id_key), ei.svlen, ej.svlen, "positions", ei.posA, ej.posA, ei.svtype,
            #          ej.svtype, 'SU=', ei.su, ej.su, ei.haplotypes, ej.haplotypes)
            #     seen.add(id_key)
            #     continue

            overlap = max(0, min(ei.posA, ej.posA) - max(ei.posB, ej.posB))
            if ei.spanning > 0 and ej.spanning > 0 and overlap == 0 and ei.svtype != "INS":
                disjoint_edges.add(id_key)
                continue

        if (same_sample and ei.svtype == "DEL" and ei.su < 3 and ej.su < 3 and
                not any_contigs_to_check and ei.spanning == 0 and ej.spanning == 0 and ei.sc == 0 and ej.sc == 0):
            continue

        ml = max(int(ei.svlen), int(ej.svlen))
        if ml == 0 and ei.svtype != 'TRA':
            continue

        l_ratio = min(int(ei.svlen), int(ej.svlen)) / ml if ml else 1

        if paired_end:  # Merge insertion-like sequences aggressively
            ei_ins_like = ei.svtype in "INSDUPTRA"
            ej_ins_like = ej.svtype in "INSDUPTRA"
            recpi_overlap = is_reciprocal_overlapping(ei.posA, ei.posB, ej.posA, ej.posB)
            if ei_ins_like and ej_ins_like and abs(ei.posA - ej.posA) < 50 and (ei.spanning == 0 and ej.spanning == 0 and ei.remap_score == 0 and ej.remap_score == 0):
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(idx, jdx, loci_same=False)
                continue
            elif recpi_overlap and ei.svtype == "DEL" and ej.svtype == "DEL" and ei.remap_score > 0 and ej.remap_score > 0:
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(idx, jdx, loci_same=False)
                continue

        # merge_conditions_met = False
        # if ins_dup or ei.svtype == "TRA" or ej.svtype == "TRA":
        #     merge_conditions_met = True
        # elif ei.svtype == "INS":
        #     if aggressive_ins_merge or paired_end: #(paired_end and isinstance(ei.variant_seq, str) and isinstance(ej.variant_seq, str) and l_ratio > 0.7):
        #         merge_conditions_met = True
        #     elif ml > 0 and l_ratio > 0.7:
        #         merge_conditions_met = True
        # else:
        #     spd = span_position_distance(ei, ej)
        #     recpi_overlap = is_reciprocal_overlapping(ei.posA, ei.posB, ej.posA, ej.posB)
        #     both_in_include = intersecter(tree, ei.chrA, ei.posA, ei.posA + 1) and intersecter(tree, ei.chrB, ei.posB, ei.posB + 1)
        #
        #     merge_conditions_met = (
        #             (paired_end and ei.spanning > 0 and ej.spanning > 0 and (recpi_overlap or spd > 0.3)) or
        #             ((recpi_overlap or spd > 0.3 or (loci_similar and any_contigs_to_check)) and not both_in_include)
        #     )


        # if not merge_conditions_met:
        #     echo("merge conditions not met", i_id, j_id, ei.svlen, ej.svlen)
        #     continue

        # if same_sample and not paired_end and not consistent_alignment_and_cigars(ei, ej, l_ratio):
        #     echo("inconsistent cigars")
        #     continue

        # echo('hi', ei.type, ej.type)

        # Loci are similar, check contig match or reciprocal overlap
        if not any_contigs_to_check:
            one_is_imprecise = (not ei.preciseA or not ei.preciseB or ei.svlen_precise or
                                not ej.preciseA or not ej.preciseB or ej.svlen_precise)
            if ml > 0 and (l_ratio > 0.5 or (one_is_imprecise and l_ratio > 0.3)):
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(idx, jdx, loci_same=False)
                # echo("MERGED2", ei.svlen, ej.svlen, ei.svtype)
                continue
            elif ei.svtype == 'TRA' or ei.svtype == 'TRA':
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(idx, jdx, loci_same=False)
        elif procs > 1:
            job.append(
                (ci, ci2, ci_alt, cj, cj2, cj_alt, ei, ej, paired_end, idx, jdx, same_sample)
            )
        else:
            edge = process_contig_aignments(ci, ci2, ci_alt, cj, cj2, cj_alt, ei, ej, paired_end, idx, jdx, same_sample)
            if not edge:
                continue
            else:
                # echo("MERGED3", ei.svlen, ej.svlen, ei.svtype)
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(edge[0], edge[1], loci_same=True)

    if job:
        for edge in pool.starmap(process_contig_aignments, job):
            if not edge:
                continue
            else:
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(edge[0], edge[1], loci_same=True)

    return G, disjoint_edges


def split_graph_by_forbidden_edges(G, forbidden_pairs, potential):

    def bfs_find_and_cut(start, target, work_graph, potential):
        queue = deque([(start, None)])  # (node, parent)
        visited = {start}
        parent_map = {}  # Keep track of how we reached each node
        while queue:
            current, parent = queue.popleft()
            for neighbor in work_graph.neighbors(current):
                if neighbor == target:
                    work_graph.remove_edge(current, neighbor)
                    return True
                if neighbor not in visited:
                    visited.add(neighbor)
                    parent_map[neighbor] = current
                    queue.append((neighbor, current))
        return False

    work_graph = G.copy()
    node_component_size = {}
    components = nx.connected_components(work_graph)

    for idx, component in enumerate(components):
        component_size = len(component)
        if len(component) > 3:  # If one node connects all others, don't cut
            max_degree = max(work_graph.degree(u) for u in component)
            if max_degree == len(component) - 1:
                continue
        for node in component:
            node_component_size[node] = component_size

    for node1, node2 in forbidden_pairs:
        if node1 not in work_graph or node2 not in work_graph:
            continue
        # For n=3, assume nodes really should be merged. This can arise when
        # one large INDEL matches two smaller INDELS that are adjacent on the same read
        # Don't cut if graph is triangle, do otherwise
        component_size = node_component_size.get(node1, 0)
        if component_size > 50:
            continue
        # echo('component size', node1, node2, component_size, work_graph.edges(node1), work_graph.edges(node2))
        if node_component_size.get(node1, 0) == 3:
            # if len(work_graph.edges(node1)) == 2 and len(work_graph.edges(node2)) == 2:
                continue
        # echo('n out edges', (node1, node2), node_component_size.get(node1, 0), len(work_graph.edges(node1)), len(work_graph.edges(node2)))
        bfs_find_and_cut(node1, node2, work_graph, potential)

    final_components = [work_graph.subgraph(c).copy() for c in nx.connected_components(work_graph)]
    return final_components


def srt_func(item1, item2):
    support1 = item1.su + (3 * item1.spanning)
    support2 = item2.su + (3 * item2.spanning)
    if item1.svtype != item2.svtype:
        # INS types for paired-end data are often incorrect compared to e.g. DUP/DEL calls, so down weight them
        if item1.svtype in 'INSTRABND' and item1.type == 'pe':
            support1 /= 2
        elif item2.svtype in 'INSTRABND' and item2.type == 'pe':
            support2 /= 2
    if support1 > support2:
        return -1
    elif support2 > support1:
        return 1

    # If supports are equal, compare svlen
    if item1.svlen > item2.svlen:
        return -1
    elif item2.svlen > item1.svlen:
        return 1
    return 0  # If everything is equal


def merge_events(potential, max_dist, tree, paired_end=False, try_rev=False, pick_best=False, add_partners=False,
                 rel_diffs=False, diffs=15, same_sample=True, debug=False, min_size=0, aggressive_ins_merge=False,
                 skip_imprecise=False, max_comparisons=100, procs=1):
    """Try and merge similar events, use overlap of both breaks points
    """
    max_dist = max_dist / 2
    if len(potential) <= 1:
        return potential
    # Cluster events on graph
    G = nx.Graph()
    G, forbidden_edges = enumerate_events(G, potential, max_dist, try_rev, tree, paired_end, rel_diffs, diffs, same_sample,
                        aggressive_ins_merge=aggressive_ins_merge,
                        debug=debug, max_comparisons=max_comparisons, procs=procs)
    # for x, y in forbidden_edges:
    #     echo('forbidden edge', (x, y), potential[x].svlen, potential[y].svlen)
    found = []
    for idx, item in enumerate(potential):  # Add singletons, non-merged
        if not G.has_node(idx):
            found.append(item)

    # Try and merge SVs with identical breaks, then merge ones with less accurate breaks - this helps prevent
    # over merging SVs that are close together
    components = split_graph_by_forbidden_edges(G, forbidden_edges, potential)
    #node_to_event = {i.event_id: i for i in potential}
    cdef int k
    for grp in components:

        best = [potential[n] for n in grp]

        best.sort(key=cmp_to_key(srt_func), reverse=False)
        w0 = best[0]

        # echo('merge groups', [i for i in best])
        # echo([(b.svlen, b.su, b.svtype) for b in best], w0.svtype, w0.svlen)
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
                if w0.qnames is not None and item.qnames is not None:
                    w0.qnames |= item.qnames
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
                        elif w0.remap_score == 0 and min_size > w0.svlen < item.svlen:
                            w0.svlen = item.svlen
                elif item.svtype == "INS" and svt in {"INS","DUP","TRA","INV"}:
                    if not spanned:
                        if item.spanning:
                            w0.svlen = item.svlen
                            w0.variant_seq = item.variant_seq
                        # elif item.svlen * 0.6 < w0.svlen < item.svlen or min_size > w0.svlen < item.svlen:
                        elif w0.remap_score == 0 and min_size > w0.svlen < item.svlen:
                            w0.svlen = item.svlen
                            w0.svtype = item.svtype
                            if best_var_seq:
                                best_var_seq = None
                                w0.variant_seq = None
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
            if not w0.partners:
                w0.partners = [i.event_id for i in best[1:]]
            else:
                w0.partners += [i.event_id for i in best[1:]]
        found.append(w0)

    return found
