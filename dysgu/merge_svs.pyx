# cython: language_level=3
from __future__ import absolute_import
import numpy as np
import random
from collections import defaultdict
import networkx as nx
from dysgu import consensus, io_funcs
from dysgu.map_set_utils cimport is_reciprocal_overlapping, EventResult
from dysgu.map_set_utils import echo
from dysgu.io_funcs import intersecter
import itertools
from cython.operator import dereference


ctypedef EventResult EventResult_t

np.random.seed(0)
random.seed(0)

def get_chrom_key(ei):
    if ei.chrA != ei.chrB:
        return min(ei.chrA, ei.chrB), max(ei.chrA, ei.chrB)
    return ei.chrA


def compare_subset(potential, int max_dist, int max_comparisons):
    tmp_list = defaultdict(list)
    cdef int idx, jdx, dist1, ci_a
    for idx in range(len(potential)):
        ei = potential[idx]
        tmp_list[get_chrom_key(ei)].append((ei.posA - ei.cipos95A - max_dist, ei.posA + ei.cipos95A + max_dist, idx))

    nc2 = {k: io_funcs.iitree(v, add_value=True) for k, v in tmp_list.items()}
    for idx in range(len(potential)):
        ei = potential[idx]
        ols = nc2[get_chrom_key(ei)].allOverlappingIntervals(ei.posA, ei.posA + 1)
        ols.remove(idx)
        if len(ols) > max_comparisons:
            random.shuffle(ols)
            ols = ols[:max_comparisons]
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


def similar_locations(intra, ei, ej):
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
            #echo("no ref overlap", (ei.contig_ref_start, ej.contig_ref_start), (ei.contig_ref_end, ej.contig_ref_end))
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
    #echo(ei.contig_cigar)
    #echo(ej.contig_cigar)
    tot_non_ref = diff_i + diff_j
    tot_sv_len = ei.svlen + ej.svlen
    ratio = tot_non_ref / tot_sv_len
    #echo("RATIO", ratio, tot_non_ref, tot_sv_len)
    if ratio > 1.8: # or ratio < 0.25:
        return False
    return True

def get_consensus_seqs(ei, ej):
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
    if ci_alt and cj_alt:
        yield ci_alt, cj_alt
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


cdef int matches_with_max_gap(str cigar, int max_gap):
    # get matching bases, and filter large alignment gaps
    cdef:
        # int idx = 0
        int matches = 0
        int num = 0  # To accumulate the number as we parse through the string
        char op
    cdef bytes t = cigar.encode('utf8')
    cdef char *c_cigar = t  # Convert Python string to char*
    # cigartuples = []
    while dereference(c_cigar):
        if b'0' <= dereference(c_cigar) <= b'9':
            # Convert digit char to int and accumulate
            num = num * 10 + (dereference(c_cigar) - 48 )
        else:
            op = dereference(c_cigar)
            if (op == b'D' or op == b'I') and num > max_gap:
                return 0
            if op == b'M':
                matches += num
            # cigartuples.append((cigar[idx], num))
            num = 0
        c_cigar += 1
        # idx += 1
    return matches

cdef bint bad_alignment(alignment, ei, ej, v):
    if not alignment:
        return True

    cdef int matches = matches_with_max_gap(alignment.cigar, 30)
    if matches == 0:
        return True

    qs = alignment.query_begin
    qe = alignment.query_end
    ts = alignment.target_begin
    te = alignment.target_end_optimal
    total_sv_length = ei.svlen + ej.svlen
    unaligned_bases = min(qs, ts) + min(len(v[0]) - qe, (len(v[1]) - te))
    if unaligned_bases > total_sv_length:
        return True


def enumerate_events(G, potential, max_dist, try_rev, tree, paired_end=False, rel_diffs=False, diffs=15,
                     same_sample=True, aggressive_ins_merge=False, debug=False, max_comparisons=20):
    event_iter = compare_all(potential) if len(potential) < 50 else compare_subset(potential, max_dist, max_comparisons)

    seen, disjoint_nodes = set(), set()
    out_edges = defaultdict(int)
    #echo("LEN potential", len(potential), paired_end)
    for ei, ej, idx, jdx in event_iter:

        i_id, j_id = ei.event_id, ej.event_id

        ins_dup = (ei.svtype == "INS" and ej.svtype == "DUP") or (ei.svtype == "DUP" and ej.svtype == "INS")

        if (i_id == j_id or (i_id, j_id) in seen or (j_id, i_id) in seen or
                out_edges[idx] > 10 or out_edges[jdx] > 10 or
                (not same_sample and ei.sample == ej.sample) or
                (ej.svtype == "DEL" and ei.svtype != "DEL") or
                (ei.svtype == "DEL" and ej.svtype != "DEL")):
                # (ei.svtype in {"INS", "DUP"} and ej.svtype not in {"INS", "DUP"}) or
                # ei.svtype != ej.svtype):

            continue

        seen.add((i_id, j_id))

        intra = ei.chrA == ej.chrA and ei.chrB == ej.chrB

        loci_similar = True
        if paired_end:
            loci_similar = similar_locations(intra, ei, ej)
            if not loci_similar:
                #echo("loci not similar")
                continue

        if not intra:
            out_edges[idx] += 1
            out_edges[jdx] += 1
            G.add_edge(i_id, j_id)
            #echo("not intra added")
            continue

        #echo("merge candidate", ei.svlen, ej.svlen, "positions", ei.posA, ej.posA, ei.svtype, ej.svtype)

        one_is_imprecise = (not ei.preciseA or not ei.preciseB or ei.svlen_precise or
                            not ej.preciseA or not ej.preciseB or ej.svlen_precise)

        any_contigs_to_check, ci, ci2, ci_alt, cj, cj2, cj_alt = get_consensus_seqs(ei, ej)

        if paired_end:
            overlap = max(0, min(ei.posA, ej.posA) - max(ei.posB, ej.posB))
            if ei.spanning > 0 and ej.spanning > 0 and overlap == 0 and ei.svtype != "INS":
                disjoint_nodes.add(i_id)
                disjoint_nodes.add(j_id)
                #echo("filetered disjoint")
                continue

        if (same_sample and ei.svtype == "DEL" and ei.su < 3 and ej.su < 3 and
                not any_contigs_to_check and ei.spanning == 0 and ej.spanning == 0 and ei.sc == 0 and ej.sc == 0):
            #echo("other1")
            continue

        ml = max(int(ei.svlen), int(ej.svlen))
        if ml == 0 and ei.svtype != 'TRA':
            #echo("other2")
            continue

        l_ratio = min(int(ei.svlen), int(ej.svlen)) / ml if ml else 1

        merge_conditions_met = False
        if ins_dup or ei.svtype == "TRA" or ej.svtype == "TRA":
            merge_conditions_met = True
        elif ei.svtype == "INS":
            if aggressive_ins_merge or paired_end: #(paired_end and isinstance(ei.variant_seq, str) and isinstance(ej.variant_seq, str) and l_ratio > 0.7):
                merge_conditions_met = True
            elif ml > 0 and l_ratio > 0.7:
                merge_conditions_met = True
        else:
            spd = span_position_distance(ei, ej)
            recpi_overlap = is_reciprocal_overlapping(ei.posA, ei.posB, ej.posA, ej.posB)
            both_in_include = intersecter(tree, ei.chrA, ei.posA, ei.posA + 1) and intersecter(tree, ei.chrB, ei.posB, ei.posB + 1)

            merge_conditions_met = (
                    (paired_end and ei.spanning > 0 and ej.spanning > 0 and (recpi_overlap or spd > 0.3)) or
                    ((recpi_overlap or spd > 0.3 or (loci_similar and any_contigs_to_check)) and not both_in_include)
            )


        if not merge_conditions_met:
            #echo("merge conditions not met", i_id, j_id, ei.svlen, ej.svlen)
            continue

        if same_sample and not paired_end and not consistent_alignment_and_cigars(ei, ej, l_ratio):
            #echo("inconsistent cigars")
            continue

        # Loci are similar, check contig match or reciprocal overlap
        if not any_contigs_to_check:
            if ml > 0 and (l_ratio > 0.5 or (one_is_imprecise and l_ratio > 0.3)):
                out_edges[idx] += 1
                out_edges[jdx] += 1
                G.add_edge(i_id, j_id, loci_same=False)
                #echo("MERGED2", ei.svlen, ej.svlen, ei.svtype)
                continue
            #echo("no contigs to check fail")
        else:
            #echo("processing contig pairs", ci, ci2, ci_alt, cj, cj2, cj_alt)
            for v in contig_pairs_iter(ci, ci2, ci_alt, cj, cj2, cj_alt):
            #    echo(v)
                if same_sample:
                    res = consensus.check_contig_match(v[0], v[1], return_int=False)
                    if bad_alignment(res, ei, ej, v):
                        continue

                    G.add_edge(i_id, j_id, loci_same=True)
                    #echo("MERGED3", ei.svlen, ej.svlen, ei.svtype)

                elif consensus.check_contig_match(v[0], v[1], return_int=True):
                    G.add_edge(i_id, j_id, loci_same=True)
                else:
                    continue
                #echo("no align", v[0], v[1])
                break





            #v = (ci_alt, cj_alt) if ci_alt and cj_alt else (ci, cj) if ci and cj else None

            # v = ((ci, cj) if (ci and cj) else (ci_alt, cj_alt)) if ci_alt and cj_alt else None
            # if v is None:
            #     echo("skipped due to no contigs")
            #     continue
            #
            # if v and same_sample:
            #
            #
            #     res = consensus.check_contig_match(v[0], v[1], return_int=False)
            #     echo("RES", res, v[0], v[1], ei.svlen, ej.svlen)
            #
            #     if res:
            #         qs, qe, ts, te, cig, _, _ = res
            #         if large_cigar_gap(cig, 30):
            #             echo("GAP too big", cig)
            #             continue
            #         total_sv_length = ei.svlen + ej.svlen
            #         aligned_bases = max(qe - qs, te - ts)
            #         unaligned_bases = qs + ts + (len(v[0]) - qe) + (len(v[1]) - te)
            #         echo(unaligned_bases, total_sv_length, aligned_bases, cig)
            #         echo(ei.contig_cigar, ej.contig_cigar)
            #         #if unaligned_bases < total_sv_length:
            #         # if aligned_bases > ref_overlap or unaligned_bases < total_sv_length:
            #         if unaligned_bases < total_sv_length:
            #             G.add_edge(i_id, j_id, loci_same=True)
            #             echo("MERGED3", ei.svlen, ej.svlen, ei.svtype)
            #             continue
            #
            # elif v and not same_sample and consensus.check_contig_match(v[0], v[1], return_int=True):
            #     G.add_edge(i_id, j_id, loci_same=True)
            #     continue
            #
            # # else:
            # # Handle contig matching for alternate sequences
            # for alt_combination in [(ci_alt, cj), (ci_alt, cj2), (cj_alt, ci), (cj_alt, ci2)]:
            #     echo("testing alt combo")
            #     if alt_combination[0] and alt_combination[1] and consensus.check_contig_match(alt_combination[0], alt_combination[1], return_int=True):
            #         out_edges[idx] += 1
            #         out_edges[jdx] += 1
            #         G.add_edge(i_id, j_id, loci_same=True)
            #         # echo(f"MERGED4")
            #         echo("MERGED4", ei.svlen, ej.svlen, ei.svtype)
            #         break

    return G, disjoint_nodes


def cut_components(G, disjoint_nodes):
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
                 skip_imprecise=False, max_comparisons=100):
    """Try and merge similar events, use overlap of both breaks points
    """
    max_dist = max_dist / 2
    if len(potential) <= 1:
        return potential
    # Cluster events on graph
    G = nx.Graph()
    G, disjoint_nodes = enumerate_events(G, potential, max_dist, try_rev, tree, paired_end, rel_diffs, diffs, same_sample,
                        aggressive_ins_merge=aggressive_ins_merge,
                        debug=debug, max_comparisons=max_comparisons)
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
        #echo(grp)

        best = [node_to_event[n] for n in grp]

        best.sort(key=srt_func, reverse=True)
        w0 = best[0]

        #echo(w0.contig, w0.contig2)
        #echo([b.svtype for b in best])
        #echo([b.svlen for b in best], w0.svlen, w0.posA, w0.svtype, w0.rep)
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
                elif item.svtype == "INS" and svt in {"INS","DUP","TRA","INV"}:
                    if not spanned:
                        if item.spanning:
                            w0.svlen = item.svlen
                            w0.variant_seq = item.variant_seq
                        # elif item.svlen * 0.6 < w0.svlen < item.svlen or min_size > w0.svlen < item.svlen:
                        elif min_size > w0.svlen < item.svlen:
                            #echo("best var seq", best_var_seq)
                            # if best_var_seq == -1:
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
            w0.partners = [i.event_id for i in best[1:]]
        found.append(w0)
    return found
