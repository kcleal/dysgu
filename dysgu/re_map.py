from skbio.alignment import StripedSmithWaterman
from dysgu.map_set_utils import is_overlapping
from dysgu.coverage import merge_intervals
from dysgu.assembler import check_contig_match, compute_rep

# from dysgu.ksw2_utils import py_ksw2_alignment

import edlib
import click


def echo(*args):
    click.echo(args, err=True)


def get_clipped_seq(cont, position, cont_ref_start, cont_ref_end):
    if cont and len(cont) > 5:
        start_i = 1
        while cont[start_i].islower():
            start_i += 1
            if start_i == len(cont):
                break
        end_i = len(cont) - 1
        while cont[end_i].islower():
            end_i -= 1
            if end_i < 0:
                break
        extra = 5
        if abs(cont_ref_start - position) < abs(cont_ref_end - position):
            if start_i > 8:
                left_clip = cont[:start_i]  # should be - 0
                return left_clip, 0, len(cont) - end_i
        else:
            if len(cont) - end_i > 8:
                right_clip = cont[end_i + 1:]  # + 1
                return right_clip, 1, start_i


def filter_bad_alignment(align, event, idx, begin, end, break_position, filter_expected=True):
    pos = event["pos" + idx]
    score = align.optimal_alignment_score
    span = align.query_end - align.query_begin + 1
    seq1 = align.aligned_query_sequence
    seq2 = align.aligned_target_sequence

    if not seq1 or not seq2:
        return -1
    if align.target_begin > 8 and len(align.target_sequence) - align.target_end_optimal > 8:
        return -1

    distance_to_break = min(abs(begin - break_position), abs(end - break_position))
    large_gap_penalty = 24
    gapped_score = score

    if distance_to_break > 200:
        gapped_score = score - large_gap_penalty

    # rep = 0
    # if distance_to_break > 50:  # Increase alignment stringency for deletions with repetitive regions
    #     rep = compute_rep(align.query_sequence)  # still penalized insertions todo fix this
    #     expect_thresh = 0.9
    # else:
    expect_thresh = 0.7

    if gapped_score > 12:
        if is_overlapping(begin - 1, end + 1, pos, pos + 1):
            return 1
        elif gapped_score > 20:
            expected = span * 2  # 2 x match score
            # if gaps at both ends of alignment increase stringency
            if align.target_begin >= 2 and align.target_end_optimal < len(align.target_sequence) - 2:
                expected = span * 4

            # if rep:
            #     expected += int(25 * rep)  # rep is float between 0 - 1
            #
            # if not filter_expected:
            #     return 1
            # echo(span, score, expected, float(score) / expected, rep)
            if span > 12 and float(score) / expected > expect_thresh:
                return expected - score
            else:
                return -1
    return -1


# def bad_soft_clip(break_position, clip_side, ref_seq_big, gstart, clip_seq):
#     # check that the soft-clip cant be aligned to the local genomic sequence using less stringent alignment params
#     if clip_side == 0:
#
#         ref_start = break_position - len(clip_seq)  # todo set as parameter
#         ref_end = break_position
#     else:
#         ref_start = break_position
#         ref_end = break_position + len(clip_seq)
#
#     start_idx = ref_start - gstart
#     start_idx = 0 if start_idx < 0 else start_idx
#     end_idx = ref_end - gstart
#
#     ref_seq_clipped = ref_seq_big[start_idx:end_idx]
#
#     if not ref_seq_clipped or ref_seq_clipped[0] in "nN" or ref_seq_clipped[-1] in "nN":
#         return True
#
#     aln = StripedSmithWaterman(ref_seq_clipped, match_score=2, mismatch_score=-3, gap_open_penalty=3, gap_extend_penalty=1)
#     a = aln(clip_seq)
#     score = a.optimal_alignment_score
#     span = a.query_end - a.query_begin + 1
#     expected = span * 2
#     if span / len(clip_seq) > 0.7 and score / expected > 0.5:
#         return True


def merge_align_regions(locations):
    # Merge any similar alignment regions found by edlib, used to get the bounds of the alignment region
    if len(locations) <= 1:
        return locations
    merge_dist = 10
    new_l = []
    for s, e in locations:
        if len(new_l) == 0:
            new_l.append([s, e])
        last = new_l[-1]

        if abs(s - last[0]) < merge_dist and abs(e - last[1]) < merge_dist:
            new_l[-1][1] = e
        else:
            return None
    return new_l


def switch_sides(e):
    # switch sides
    chrA, posA, cipos95A, contig2  = e["chrA"], e["posA"], e["cipos95A"], e["contig2"]
    contig2_ref_start, contig2_ref_end, contig2_left_weight, contig2_right_weight = e["contig_ref_start"], e["contig_ref_end"], e["contig_left_weight"], e["contig_right_weight"]
    # e["chrA"] = e["chrB"] assumed to be the same
    e["posA"] = e["posB"]
    e["cipos95A"] = e["cipos95B"]
    e["chrB"] = chrA
    e["posB"] = posA
    e["cipos95B"] = cipos95A

    e["contig2"] = e["contig"]
    e["contig"] = contig2

    e["contig_ref_start"] = e["contig2_ref_start"]
    e["contig2_ref_start"] = contig2_ref_start

    e["contig_ref_end"] = e["contig2_ref_end"]
    e["contig2_ref_end"] = contig2_ref_end

    e["contig_left_weight"] = e["contig2_left_weight"]
    e["contig2_left_weight"] = contig2_left_weight

    e["contig_right_weight"] = e["contig2_right_weight"]
    e["contig2_right_weight"] = contig2_right_weight
    return e


def remap_soft_clips(events, ref_genome, min_sv_len, input_bam, keep_unmapped=True):

    new_events = []
    ref_locs = []

    for count, e in enumerate(events):

        e["remapped"] = 0
        e["remap_score"] = 0
        e["remap_ed"] = 0
        if 'svlen_precise' not in e:
            e['svlen_precise'] = 1

        if e["chrA"] != e["chrB"]:
            new_events.append(e)
            continue

        try_remap = False
        if (e["contig"] or e["contig2"]) and (e["svlen"] < 1000):
            if not e['svlen_precise']:
                try_remap = True

        if not try_remap:
            e["modified"] = 0
            new_events.append(e)
            continue

        else:
            if e["posA"] <= e["posB"]:
                ref_locs.append((e["chrA"], e["posA"], e["posB"], count))
            else:
                ref_locs.append((e["chrA"], e["posB"], e["posA"], count))

    for chrom, gstart, gend, grp_idxs in merge_intervals(ref_locs, pad=1500, add_indexes=True):
        if gstart < 0:
            gstart = 0
        try:
            ref_seq_big = ref_genome.fetch(chrom, gstart, gend).upper()
        except (ValueError, KeyError, IndexError) as errors:
            # Might be different reference genome version, compared to bam genome
            continue

        # note this doesnt parellelize well with multiprocessing.pool, suspect serializing is too slow
        # could try and parallel before the fetch command above using publish/subscribe model
        for index in grp_idxs:

            e = events[index]
            added = 0
            high_quality_clip = False
            skip_event = False
            for cont, idx in (("contig", "A"), ("contig2", "B")):
                if cont in e and e[cont]:

                    break_position = e["pos" + idx]
                    clip_res = get_clipped_seq(e[cont], break_position, e[cont + "_ref_start"], e[cont + "_ref_end"])

                    if not clip_res:
                        continue
                    clip_seq, clip_side, length_other_clip = clip_res

                    if length_other_clip > 3 and e['ref_bases'] < 50:
                        continue

                    if clip_side == 0:
                        w = e[cont + "_left_weight"]
                        if not w > 10:  # todo set as a parameter option
                            continue
                        elif w > 400:
                            high_quality_clip = True
                    else:
                        w = e[cont + "_right_weight"]
                        if not w > 10:
                            continue
                        elif w > 400:
                            high_quality_clip = True

                    # if bad_soft_clip(break_position, clip_side, ref_seq_big, gstart, clip_seq):
                    #     continue

                    ref_start = break_position - 500  # todo set as parameter
                    ref_end = break_position + 500

                    start_idx = ref_start - gstart
                    start_idx = 0 if start_idx < 0 else start_idx
                    end_idx = ref_end - gstart

                    ref_seq_clipped = ref_seq_big[start_idx:end_idx]

                    ref_seq_start = gstart + start_idx
                    # ref_seq_end = gstart + end_idx
                    if not ref_seq_clipped or ref_seq_clipped[0] in "nN" or ref_seq_clipped[-1] in "nN":
                        skip_event = True
                        break

                    # Large alignment region

                    el = edlib.align(clip_seq.upper(), ref_seq_clipped, mode="HW", task="locations")
                    locs = merge_align_regions(el['locations'])

                    if not locs:
                        continue

                    l_start, l_end = locs[0]

                    ref_start2 = ref_seq_start + l_start

                    # if clip_side == 1 and e[cont + "_ref_start"] < ref_start2:
                    #     l_start -= (ref_start2 - e[cont + "_ref_start"])
                    #     echo(ref_start2 - e[cont + "_ref_start"])

                    ref_seq2 = ref_seq_clipped[l_start:l_end+1]

                    aln = StripedSmithWaterman(ref_seq2, match_score=2, mismatch_score=-8, gap_open_penalty=12, gap_extend_penalty=1)
                    a = aln(clip_seq)
                    # echo(clip_seq)
                    # echo(a)
                    # echo(locs)
                    # echo(el)
                    # a = py_ksw2_alignment(ref_seq2, clip_seq, 2, -8, 12, 1)
                    # echo(a)
                    # echo(ref_start2)
                    # echo(e["posA"], e["posB"])
                    # echo(ref_seq2)
                    # echo(e["contig_ref_start"])
                    # echo(clip_side)

                    score = a.optimal_alignment_score

                    aln_q_end = a.query_end
                    aln_q_begin = a.query_begin
                    aln_t_begin = a.target_begin
                    target_end_optimal = a.target_end_optimal

                    q_begin = ref_start2 + aln_q_begin
                    q_end = ref_start2 + aln_q_end

                    edit_dist = filter_bad_alignment(a, e, idx, q_begin, q_end, break_position, filter_expected=True)
                    # echo("filter passed", edit_dist != -1)
                    # if edit_dist == -1:
                    #
                    #     a = py_ksw2_alignment(ref_seq2, clip_seq, 2, -8, 12, 1)
                    #     # aln = StripedSmithWaterman(ref_seq2, match_score=2, mismatch_score=-8, gap_open_penalty=12,
                    #     #                            gap_extend_penalty=1)
                    #     # a = aln(clip_seq)
                    #
                    #     score = a.optimal_alignment_score
                    #
                    #     aln_q_end = a.query_end
                    #     aln_q_begin = a.query_begin
                    #     aln_t_begin = a.target_begin
                    #     target_end_optimal = a.target_end_optimal
                    #
                    #     q_begin = ref_start2 + aln_q_begin
                    #     q_end = ref_start2 + aln_q_end
                    #
                    #     edit_dist = filter_bad_alignment(a, e, idx, q_begin, q_end, break_position)

                        # if edit_dist != -1:
                        #     echo("2", e)
                        #     echo(a.cigar)
                        #     echo(a)
                        #     echo()

                    if edit_dist != -1:

                        pos = e["pos" + idx]
                        if clip_side == 0:
                            if q_end + 1 >= pos:
                                kind = "INS"
                                break_point = pos
                                break_point2 = pos
                                overlap = q_end - pos
                                svlen = len(clip_seq) - target_end_optimal + overlap
                            else:
                                ref_gap = pos - q_end
                                target_gap = len(clip_seq) - target_end_optimal

                                if target_gap > ref_gap:
                                    kind = "INS"
                                    break_point = pos
                                    break_point2 = pos
                                    svlen = target_gap

                                else:
                                    kind = "DEL"
                                    break_point = pos
                                    break_point2 = q_end
                                    svlen = ref_gap

                            # discard alignments with large unmapped overhang
                            if aln_t_begin > svlen:
                                continue
                        else:

                            if q_begin - 1 <= pos:
                                kind = "INS"
                                break_point = pos
                                break_point2 = pos
                                if q_end > pos:
                                    svlen = pos - q_begin + aln_t_begin
                                else:
                                    svlen = max(q_end, pos) -  min(q_begin, pos)

                            else:
                                ref_gap = q_begin - pos
                                target_gap = aln_t_begin
                                if target_gap > ref_gap:
                                    kind = "INS"
                                    break_point = pos
                                    break_point2 = pos
                                    svlen = target_gap
                                else:
                                    kind = "DEL"
                                    break_point = pos
                                    break_point2 = q_begin
                                    svlen = abs(break_point2 - break_point)

                            if len(clip_seq) - target_end_optimal > svlen:
                                continue

                        if svlen < min_sv_len:
                            continue

                        if kind == "DEL":
                            span = a.query_end - a.query_begin + 1
                            if span < len(clip_seq) * 0.4 and span < 50:
                                continue

                        if abs(svlen - e['svlen']) > 20:
                            e["remap_ed"] = edit_dist
                            e["remapped"] = 1
                            e["remap_score"] = score
                            e['svtype'] = kind
                            e['svlen'] = svlen
                            e['pos' + idx] = break_point
                            if idx == "A":
                                other = "B"
                            else:
                                other = "A"
                            e['pos' + other] = break_point2
                            e['cipos95A'] = 0
                            e['cipos95B'] = 0

                            # switch if nessasary
                            if e['posA'] > e['posB']:
                                e = switch_sides(e)

                            new_events.append(e)
                            added = 1
                            break  # dont analyse contig2

                if added:
                    break

            if not added and not skip_event and high_quality_clip and keep_unmapped:
                new_events.append(e)

    return new_events
    # check for smaller indels that could explain the SV
    # checked = []

    # for e in new_events:
    #     passed = True
    #     if e["svtype"] == "DEL":
    #         # get local reads with biggest dels and check if found del can be aligned
    #         reads_with_dels = []
    #         for a in input_bam.fetch(e["chrA"], e["posA"], e["posA"] + 1):
    #             if a.cigartuples is None:
    #                 continue
    #             ds = sorted([i for i in a.cigartuples if i[0] == 2 and i[1] < 30], key=lambda x: x[1], reverse=True)
    #             if len(ds):
    #                 if ds[0][1] > 10:
    #                     passed = False
    #                     break
            #         reads_with_dels.append((ds[0], a.seq))
            # if len(reads_with_dels) == 0:
            #     continue
            # srt_rd = sorted(reads_with_dels, reverse=True)
            #
            # clip_seq = None
            # if e["contig"]:
            #     clip_res = get_clipped_seq(e[cont], e["posA"], e["contig_ref_start"], e["contig_ref_end"])
            #     if clip_res:
            #         clip_seq, _, _ = clip_res
            # else:
            #     clip_res = get_clipped_seq(e[cont], e["posB"], e["contig2_ref_start"], e["contig2_ref_end"])
            #     if clip_res:
            #         clip_seq, _, _ = clip_res
            #
            # for _, seq in srt_rd[:3]:
            #
            #     if check_contig_match(seq, clip_seq, return_int=True):
            #         passed = False
            #         break
    #     if passed:
    #         checked.append(e)
    #
    # return checked