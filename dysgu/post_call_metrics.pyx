#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

from collections import defaultdict
import numpy as np
cimport numpy as np
from pysam.libchtslib cimport bam_get_qname

from dysgu.map_set_utils cimport unordered_map
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement

from libc.math cimport fabs as c_fabs
from libc.stdint cimport uint64_t
from libcpp.vector cimport vector as cpp_vector

from dysgu.map_set_utils import echo

import math


cdef float soft_clip_qual_corr(reads):
    """Function to compute the correlation of quality values between the soft-clipped portions of reads. values of
    1 imply no correlation whereas closer to 0 implies that quality values are correlated"""
    cdef const unsigned char[:] quals
    cdef int idx, x, i, j
    cdef unordered_map[int, cpp_vector[int]] qq
    for r in reads:

        quals = r.query_qualities
        if r.cigartuples[0][0] == 4:
            idx = 0
            for x in range(r.pos - r.cigartuples[0][1], r.pos):
                qq[x].push_back(quals[idx])
                idx += 1
        if r.cigartuples[-1][0] == 4:
            end_pos = r.query_alignment_end + r.cigartuples[-1][1]
            for idx in range(-1, - r.cigartuples[-1][1], -1):
                qq[end_pos].push_back(quals[idx])
                end_pos -= 1

    if qq.empty():
        return -1

    cdef float all_z = 0
    cdef bint seen = False
    cdef cpp_vector[int] all_v
    cdef float sum_all_v = 0
    cdef cpp_vector[int] second
    cdef int first
    cdef float mean_second, sum_second, sum_diff
    cdef unordered_map[int, cpp_vector[int]].iterator qq_iter

    with nogil:
        qq_iter = qq.begin()
        while qq_iter != qq.end():

            second = dereference(qq_iter).second

            if second.size() > 1:
                sum_second = 0
                for i in range(second.size()):
                    sum_second += second[i]
                mean_second = sum_second / float(second.size())
                for i in range(second.size()):
                    all_z += c_fabs(second[i] - mean_second)

                sum_all_v += sum_second
                all_v.insert(all_v.end(), second.begin(), second.end())  # concat

                seen = True
            preincrement(qq_iter)

    if not seen:
        return -1

    cdef float mean_all_v = sum_all_v / float(all_v.size())

    cdef float z = 0
    for i in range(all_v.size()):
        z += c_fabs(all_v[i] - mean_all_v)

    if z == 0:
        return -1

    return all_z / z


# stolen from svtyper with minor modifications
# https://github.com/hall-lab/svtyper/blob/master/svtyper/singlesample.py
# efficient combinatorial function to handle extremely large numbers
def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in range(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# return the genotype and log10 p-value
def bayes_gt(ref, alt, is_dup):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [1e-2, 0.2, 1/3.0]
    else:
        p_alt = [1e-3, 0.5, 0.9]

    total = ref + alt
    log_combo = log_choose(total, alt)

    lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return lp_homref, lp_het, lp_homalt


def get_gt_metric(events, infile):

    new_events = []
    pad = 5
    for e in events:

        if infile is None:
            e["GT"] = "./."
            e["GQ"] = -1
            continue

        # d = set([])
        # for r in infile.fetch(e["chrA"], e["posA"] - pad, e["posA"] + pad):
        #     unique_name = (r.qname, r.flag, r.pos, r.rname)
        #     d.add(unique_name)
        # if e["chrA"] != e["chrB"] or (e["posB"] != e["posA"] and abs(e["posB"] - e["posA"]) > pad*2):
        #     for r in infile.fetch(e["chrA"], e["posA"] - pad, e["posA"] + pad):
        #         unique_name = (r.qname, r.flag, r.pos, r.rname)
        #         d.add(unique_name)

        # total = len(d)
        support_reads = e["su"] - e["spanning"]  # spanning are counted twice
        ref = int(e["raw_reads_10kb"]) #total - support_reads
        if ref < 0:
            ref = 0
        # echo(ref, support_reads)
        gt_lplist = bayes_gt(ref, support_reads, e["svtype"] == "DUP")
        best, second_best = sorted([ (i, e) for i, e in enumerate(gt_lplist) ], key=lambda x: x[1], reverse=True)[0:2]
        gt_idx = best[0]

        gt_sum = 0
        result = {}
        for gt in gt_lplist:
            try:
                gt_sum += 10**gt
            except OverflowError:
                gt_sum += 0
        if gt_sum > 0:
            gt_sum_log = math.log(gt_sum, 10)
            sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
            phred_gq = min(-10 * (second_best[1] - best[1]), 200)
            result['GQ'] = int(phred_gq)
            result['SQ'] = sample_qual
            if gt_idx == 1:
                result['GT'] = '0/1'
            elif gt_idx == 2:
                result['GT'] = '1/1'
            elif gt_idx == 0:
                result['GT'] = '0/0'
        else:
            result['GQ'] = '.'
            result['SQ'] = '.'
            result['GT'] = './.'
        e.update(result)
        # echo(result)

    return events

    # return new_events