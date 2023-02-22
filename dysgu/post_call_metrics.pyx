#cython: language_level=3, boundscheck=True, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

import logging
import numpy as np
cimport numpy as np

from dysgu.map_set_utils cimport unordered_map, EventResult
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement

from libc.math cimport fabs as c_fabs
from libcpp.vector cimport vector as cpp_vector

from dysgu.map_set_utils import echo
from dysgu import re_map
from dysgu.io_funcs import reverse_complement, intersecter
from dysgu.assembler import compute_rep

import zlib
import math
import array
import pickle
import glob
import gzip
import pandas as pd
import warnings

pd.options.mode.chained_assignment = None

from dysgu.scikitbio._ssw_wrapper import StripedSmithWaterman
import time
from libc.stdint cimport int16_t
from sys import byteorder
import os

ctypedef EventResult EventResult_t


class BadClipCounter:

    def __init__(self, n_references, low_mem, temp_dir):
        self.low_mem = low_mem
        self.temp_dir = temp_dir
        self.n_refs = n_references
        if not self.low_mem:
            self.clip_pos_arr = [array.array("L", []) for i in range(n_references)]  # 32 bit unsigned int
        else:
            self.clip_pos_arr = [open(f"{self.temp_dir}/{i}.badclip.bin", "wb") for i in range(n_references)]

    def tidy(self):
        if self.low_mem:
            [os.remove(f"{self.temp_dir}/{i}.badclip.bin") for i in range(self.n_refs)]

    def add(self, chrom, position):
        if not self.low_mem:
            self.clip_pos_arr[chrom].append(position)
        else:
            self.clip_pos_arr[chrom].write(position.to_bytes(4, byteorder))

    def sort_arrays(self):

        if not self.low_mem:
            self.clip_pos_arr = [np.array(i, dtype=np.dtype("i")) for i in self.clip_pos_arr]
        else:
            [i.close() for i in self.clip_pos_arr]
            self.clip_pos_arr = [np.fromfile(f"{self.temp_dir}/{i}.badclip.bin", np.int32) for i in range(self.n_refs)]

        [i.sort() for i in self.clip_pos_arr]

    def count_near(self, int chrom, int start, int end):

        py_a = self.clip_pos_arr[chrom]
        cdef int len_a = len(py_a)
        if chrom < 0 or chrom >= len_a or len_a == 0:
            return 0

        cdef int[:] a = py_a

        if start < 0:
            start = 0

        cdef int i = np.searchsorted(py_a, start)
        if i < 0:
            i = 0
        if i >= len_a - 1:
            return 0

        # search forward and backwards
        cdef int count = 0
        cdef int p

        while True:
            if i >= len_a:
                break
            p = a[i]
            if p <= end:
                count += 1
                i += 1
            else:
                break

        return count


def get_badclip_metric(events, bad_clip_counter, bam, regions):

    bad_clip_counter.sort_arrays()
    new_events = []
    cdef EventResult_t e
    for e in events:
        count = 0
        if e.chrA == e.chrB and abs(e.posB - e.posA) < 500:

            start = min(int(e.posA), int(e.posB)) - 500
            end = max(int(e.posA), int(e.posB)) + 500
            if not intersecter(regions, e.chrA, e.posA, e.posA + 1) and \
                not intersecter(regions, e.chrB, e.posB, e.posB + 1):

                count = bad_clip_counter.count_near(bam.gettid(e.chrA), start, end)
        else:
            c1, c2 = 0, 0
            if not intersecter(regions, e.chrA, e.posA, e.posA + 1):
                c1 = bad_clip_counter.count_near(bam.gettid(e.chrA), e.posA - 500, e.posA + 500)
            if not intersecter(regions, e.chrB, e.posB, e.posB + 1):
                c2 = bad_clip_counter.count_near(bam.gettid(e.chrB), e.posB - 500, e.posB + 500)
            count = c1 + c2

        e.bad_clip_count = count
        e.ras = 0
        e.fas = 0

        if e.spanning > 0:
            new_events.append(e)
            continue

        clip_res = None

        max_score_rev = 0
        max_score_forward = 0

        if e.contig and len(e.contig) < 1000:
            break_position = e.posA
            clip_res = re_map.get_clipped_seq(e.contig, e.posA, e.contig_ref_start, e.contig_ref_end)
            if clip_res:
                fc = clip_res[0]
                rc = reverse_complement(fc, len(fc))
                query = StripedSmithWaterman(e.contig, gap_extend_penalty=1, suppress_sequences=True)
                align = query(rc)

                if align.optimal_alignment_score > max_score_rev:
                    max_score_rev = align.optimal_alignment_score

                align = query(fc)
                if align.optimal_alignment_score > max_score_forward:
                    max_score_forward = align.optimal_alignment_score

        if e.contig2 and len(e.contig2) < 1000:
            break_position = e.posB
            clip_res = re_map.get_clipped_seq(e.contig2, e.posB, e.contig2_ref_start, e.contig2_ref_end)
            if clip_res:
                fc = clip_res[0]
                rc = reverse_complement(fc, len(fc))
                query = StripedSmithWaterman(e.contig2, gap_extend_penalty=1, suppress_sequences=True)
                align = query(rc)

                if align.optimal_alignment_score > max_score_rev:
                    max_score_rev = align.optimal_alignment_score

                align = query(fc)
                if align.optimal_alignment_score > max_score_forward:
                    max_score_forward = align.optimal_alignment_score

        e.ras = max_score_rev
        e.fas = max_score_forward

        new_events.append(e)

    return new_events


cdef float soft_clip_qual_corr(reads):
    """Function to compute the correlation of quality values between the soft-clipped portions of reads. values of
    1 imply no correlation whereas closer to 0 implies that quality values are correlated"""
    cdef const unsigned char[:] quals
    cdef int idx, x, i, j
    cdef unordered_map[int, cpp_vector[int]] qq

    for r in reads:

        quals = r.query_qualities
        if r.cigartuples is None or quals is None:
            continue
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


class CoverageAnalyser(object):

    def __init__(self, temp_dir):

        self.pad_size = 1000
        self.temp_dir = temp_dir
        self.chrom_cov_arrays = {}

        chrom_cov_arrays = {}
        if os.path.exists(self.temp_dir):
            for pth in glob.glob(self.temp_dir + "/*.dysgu_chrom.bin"):
                chrom_name = pth.split("/")[-1].split(".")[0]
                chrom_cov_arrays[chrom_name] = np.fromfile(pth, dtype="int16")
            if chrom_cov_arrays:
                self.chrom_cov_arrays = chrom_cov_arrays
                logging.info("Loaded n={} chromosome coverage arrays from {}".format(len(chrom_cov_arrays), self.temp_dir))
        else:
            logging.warning("Coverage track not loaded, working directory does not exist {}".format(self.temp_dir))

    def process_events(self, events):

        # Try and load coverage tracks
        cdef EventResult_t e
        for e in events:
            e.fcc = -1
            e.inner_cn = -1
            e.outer_cn = -1

        for e in events:
            t0 = time.time()
            if e.svtype == "DEL" or e.svtype == "DUP" or e.svtype == "INV":
                if e.svlen > 2000 or abs(e.posB - e.posA) > 2000:
                    self.process_two_windows(e)
                else:
                    self.process_one_window(e)
            elif e.svtype == "INS":
                self.process_insertion(e)
            else:  # TRA
                self.process_two_windows(e)

        return events

    def process_one_window(self, EventResult_t e):
        cdef int start, end, s
        cdef float left, right, middle, sides, fc
        start, end = e.posA, e.posB
        if end < start:
            s = start
            start = end
            end = s
        if self.chrom_cov_arrays:
            if e.chrA in self.chrom_cov_arrays:
                arr = self.chrom_cov_arrays[e.chrA]
                starting_pos = 0
            else:
                return -1
        else:
            return -1

        left = median(arr, start - starting_pos - self.pad_size, start - starting_pos)
        right = median(arr, end - starting_pos, end - starting_pos + self.pad_size)
        middle = median(arr, start - starting_pos, end - starting_pos)
        if left == -1 or right == -1 or middle == -1:
            return -1

        sides = (left + right) / 2
        fc = 0
        if e.svtype == "DEL":
            if sides == 0:
                fc = 0
            else:
                fc = middle / sides
            e.outer_cn = sides
            e.inner_cn = middle

        else:  # DUP, INV
            if middle == 0:
                fc = 0
            else:
                fc = sides / middle
            e.outer_cn = middle
            e.inner_cn = sides

        e.fcc = fc

    def process_two_windows(self, EventResult_t e):
        cdef int start, end, s
        cdef float left1, right1, left2, right2, middle, sides, fc
        start, end = e.posA, e.posB
        if e.chrA == e.chrB and end < start:
            s = start
            start = end
            end = s
        arr1, arr2 = None, None
        if self.chrom_cov_arrays:
            try:
                arr1 = self.chrom_cov_arrays[e.chrA]
                if e.chrA == e.chrB:
                    arr2 = arr1
                else:
                    arr2 = self.chrom_cov_arrays[e.chrB]
            except KeyError:
                return -1
        else:
            return -1

        left1 = median(arr1, start - self.pad_size, start)
        right1 = median(arr1, start, start + self.pad_size)

        left2 = median(arr2, end - self.pad_size, end)
        right2 = median(arr2, end, end + self.pad_size)

        if left1 == -1 or left2 == -1 or right1 == -1 or right2 == -1:
            return -1

        middle = (right1 + left2) / 2.
        sides = (left1 + right2) / 2.

        fc = 0
        if e.svtype == "DEL":
            if sides == 0:
                fc = 0
            else:
                fc = middle / sides
            e.outer_cn = sides
            e.inner_cn = middle

        elif e.svtype == "DUP":  # DUP, INV
            if middle == 0:
                fc = 0
            else:
                fc = sides / middle
            e.outer_cn = middle
            e.inner_cn = sides

        else:  # TRA
            if e.join_type == "3to5":
                e.outer_cn = sides
                e.inner_cn = middle
            elif e.join_type == "5to3":
                e.outer_cn = middle
                e.inner_cn = sides
            elif e.join_type == "3to3":
                e.outer_cn = ((left1 + left2) / 2.)
                e.inner_cn = ((right1 + right2) / 2.)
            else:  # 5to5
                e.outer_cn = ((right1 + right2) / 2.)
                e.inner_cn = ((left1 + left2) / 2.)

        e.fcc = fc

    def process_insertion(self, EventResult_t e):
        cdef float left, left_svlen, right, right_svlen, fcc, inner, outer
        if e.svlen > 10000:
            return -1

        if self.chrom_cov_arrays:
            if e.chrA in self.chrom_cov_arrays:
                arr = self.chrom_cov_arrays[e.chrA]
                starting_pos = 0
            else:
                return -1
        else:
            return -1

        pad = e.svlen if e.svlen_precise else 100

        left = median(arr, e.posA - starting_pos - self.pad_size, e.posA - starting_pos)
        left_svlen = median(arr, e.posA - pad - starting_pos, e.posA - starting_pos)

        right = median(arr, e.posA - starting_pos, e.posA - starting_pos + self.pad_size)
        right_svlen = median(arr, e.posA - starting_pos, e.posA - starting_pos + pad)

        fcc = -1
        inner = 1
        outer = 1
        if left_svlen > right_svlen and left_svlen > 0 and right != -1:
            fcc = right / left_svlen
            inner = left_svlen
            outer = right
        elif right_svlen > 0 and left != -1:
            fcc = left / right_svlen
            inner = right_svlen
            outer = left

        e.outer_cn = outer
        e.inner_cn = inner
        e.fcc = fcc

    def _get_cov(self, cn, chrom_a, chrom_b, chrom_medians):
        cdef float m, m1, m2
        if chrom_a == chrom_b:
            if chrom_medians[chrom_a]:
                m = cn / chrom_medians[chrom_a]
            else:
                m = -1
                logging.warning("Chromosome median {}: {}".format(chrom_a, chrom_medians[chrom_a]))
        else:
            m1 = chrom_medians[chrom_a]
            m2 = chrom_medians[chrom_b]
            if m1 and m2:
                m = cn / ((m1 + m2) / 2)
            else:
                m = 0
                logging.warning("Chromosome median {}: {}, {}: {}".format(chrom_a, chrom_medians[chrom_a],
                                                                          chrom_b, chrom_medians[chrom_b]))
        return m

    def normalize_coverage_values(self, events):

        if not self.chrom_cov_arrays:
            return events

        chrom_medians = {k: np.median(v[v > 0]) for k, v in self.chrom_cov_arrays.items()}
        cdef EventResult_t e
        for e in events:
            if e.outer_cn > 0:
                e.outer_cn = self._get_cov(e.outer_cn, e.chrA, e.chrB, chrom_medians)
            if e.inner_cn > 0:
                e.inner_cn = self._get_cov(e.inner_cn, e.chrA, e.chrB, chrom_medians)

        return events


def ref_repetitiveness(events, mode, ref_genome):
    cdef EventResult_t e
    for e in events:
        e.ref_rep = 0
        if e.svlen < 150 and e.svtype == "DEL":
            try:
                ref_seq = ref_genome.fetch(e.chrA, e.posA, e.posB).upper()
            except ValueError:  # todo out or range, needs fixing
                continue
            e.ref_rep = compute_rep(ref_seq)

    return events


cdef float median(np.ndarray[int16_t, ndim=1]  arr, int start, int end):
    s = int(start / 10)
    e = int(end / 10)
    if s > e:
        s1 = s
        s = e
        e = s1
    if s < 1:
        s = 1
    if e > len(arr):
        e = len(arr)
    if e > s:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            m = np.median(arr[s:e])
            if m == m:
                return m
    return -1


def nCk(n,k):
    f = math.factorial
    return f(n) // f(k) // f(n-k)


def binom_prob(n, k, p):
    # If n is > 1000, set it to 1000 so our Factorial doesn't blow up
    n = n if n <= 1000 else 1000
    binom = 0
    q = 1-p
    while k < n:
        binom += nCk(n, k) * q**(n-k) * p**k
        k += 1
    return binom


def strand_binom_t(events):
    # perform a binomial test for strand bias
    cdef EventResult_t e
    cdef int n, k
    for e in events:
        k = max(e.plus, e.minus)
        if k > 0:
            n = e.plus + e.minus
            e.strand_binom_t = binom_prob(n, k, 0.5)
        else:
            e.strand_binom_t = 1
    return events


def get_ref_base(events, ref_genome):
    cdef EventResult_t e
    for e in events:
        if not e.ref_seq:
            if e.posA == 0:
                e.posA = 1
            try:
                base = ref_genome.fetch(e.chrA, e.posA - 1, e.posA).upper()
                e.ref_seq = base
            except:
                pass

    return events


# from svtyper with minor modifications
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


def get_gt_metric(events, add_gt=False):
    cdef EventResult_t e
    if add_gt:
        logging.info("Adding genotype")
    else:
        for e in events:
            e.GQ = '.'
            # e.SQ = '.'
            e.GT = './.'
        return events


    new_events = []
    pad = 50
    for e in events:

        sup = e.su - e.spanning  # spanning are counted twice
        ocn = e.outer_cn
        icn = e.inner_cn

        if any((ocn < 0, icn < 0, sup < 0)):
            e.GQ = '.'
            # e.SQ = '.'
            e.GT = './.'
            continue

        higher_cn = max(ocn, icn)
        lower_cn = min(icn, ocn)

        support_reads = max(int((higher_cn - lower_cn)), sup)
        ref = int(lower_cn)
        if e.svtype == "INS":
            ref = int(higher_cn - sup)
            if ref < 0:
                ref = 0

        gt_lplist = bayes_gt(ref, support_reads, e.svtype == "DUP")
        best, second_best = sorted([(i, j) for i, j in enumerate(gt_lplist)], key=lambda x: x[1], reverse=True)[0:2]
        gt_idx = best[0]

        gt_sum = 0

        for gt in gt_lplist:
            try:
                gt_sum += 10**gt
            except OverflowError:
                gt_sum += 0
        if gt_sum > 0:
            gt_sum_log = math.log(gt_sum, 10)
            sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
            phred_gq = min(-10 * (second_best[1] - best[1]), 200)
            e.GQ = int(phred_gq)
            # e.SQ = sample_qual
            if gt_idx == 1:
                e.GT = '0/1'
            elif gt_idx == 2:
                e.GT = '1/1'
            elif gt_idx == 0:
                e.GT = '0/0'
        else:
            e.GQ = '.'
            # e.SQ = '.'
            e.GT = './.'

    return events


cdef del_like(EventResult_t r):

    ocn = r.outer_cn
    icn = r.inner_cn
    lower_cn = min(ocn, icn)
    t = None
    if r.inner_cn > r.outer_cn * 1.25:
        try:
            t = ((r.inner_cn - r.outer_cn) / (r.inner_cn + r.outer_cn)) * r.outer_cn
        except ZeroDivisionError:
            t = None

    if bool(r.contig) ^ bool(r.contig2):
        m = 2. * (r.sc + r.bnd - r.supp)
    else:
        m = r.sc + r.bnd - r.supp
    mr = lower_cn - m
    mr = 0 if mr < 0 else mr
    variant_supp = 0
    ref_supp = 0
    if r.pe > 0 or r.spanning > 0:
        variant_supp = r.pe + r.spanning * 0.5
        ref_supp = lower_cn

    elif r.remap_score > 0:
        if r.bnd == 0:
            variant_supp = max(r.sc, r.supp)
        else:
            if bool(r.contig) ^ bool(r.contig2):
                m = 2.5 * (r.sc + r.bnd - r.supp)
            else:
                m = 3.5 * (r.sc + r.bnd - r.supp)
            variant_supp = m
        ref_supp = lower_cn

    elif r.supp > 0:
        variant_supp = r.supp
        if t is not None:
            ref_supp = t
        else:
            ref_supp = mr

    elif r.NP > 0:
        variant_supp = r.NP * 2
        ref_supp = lower_cn
    elif r.inner_cn != -1 and r.inner_cn < r.su:
        variant_supp = r.su
        ref_supp = r.inner_cn
    elif r.outer_cn != -1 and r.inner_cn != -1:
        variant_supp = r.inner_cn
        ref_supp = r.outer_cn
    else:
        return 0, 0
        # logging.warning("No genotype for", r.chrA, r.posA, "pe", r.pe, "np", r.NP, "su", r.su, "sr", r.supp, "bnd", r.bnd, "rms", r.remap_score, "icn", r.inner_cn, "ocn", r.outer_cn,
        #       r.svtype)

    return round(ref_supp), round(variant_supp)


cdef ins_like(EventResult_t r):
    ocn = r.outer_cn
    icn = r.inner_cn
    lower_cn = min(ocn, icn)

    t = None
    if r.outer_cn > 0:
        cn_ratio = r.inner_cn / r.outer_cn
    else:
        cn_ratio = 2

    if cn_ratio > 1.25:
        try:
            t = ((r.inner_cn - r.outer_cn) / (r.inner_cn + r.outer_cn)) * r.outer_cn
        except ZeroDivisionError:
            t = None

    if bool(r.contig) ^ bool(r.contig2):
        m = 2. * (r.sc + r.bnd - r.supp)
    else:
        m = r.sc + r.bnd - r.supp
    mr = lower_cn - m
    mr = 0 if mr < 0 else mr

    variant_support = 0
    ref_support = 0

    if cn_ratio > 1.75:
        ref_support = max(0, r.inner_cn - (2. * (r.inner_cn - r.outer_cn)))
        variant_support = r.su

    elif r.pe > 0 or r.spanning > 0:
        variant_supp = r.pe + r.spanning * 0.5
        ref_support = lower_cn

    elif r.remap_score > 0:
        variant_support = m
        ref_support = mr

    elif r.supp > 0:
        variant_support = r.supp
        if t is not None:
            ref_support = t
        else:
            ref_support = mr

    elif r.NP > 0:
        variant_support = r.NP * 2
        ref_support = lower_cn - variant_support

    elif r.inner_cn > r.su:
        variant_support = r.su
        ref_support = r.inner_cn - r.su
    else:
        return 0, 0
        # logging.warning("No genotype for", r.chrA, r.posA, "pe", r.pe, "np", r.NP, "su", r.su, "sr", r.supp, "bnd", r.bnd, "rms", r.remap_score, "icn", r.inner_cn, "ocn", r.outer_cn,
        #       r.svtype)

    return round(ref_support), round(variant_support)


cdef ins_like_non_pe(EventResult_t r):
    ocn = r.outer_cn
    icn = r.inner_cn
    lower_cn = min(ocn, icn)
    higher_cn = max(ocn, icn)
    sup = r.su - r.spanning
    support_reads = max(int((higher_cn - lower_cn)), sup)
    ref = int(higher_cn - sup)
    if ref < 0:
        ref = 0
    return round(ref), round(support_reads)


def bayes_gt2(ref, alt, is_dup, a, b):
    if is_dup:  # specialized logic to handle non-destructive events such as duplications
        p_alt = b
    else:
        p_alt = a

    total = ref + alt
    log_combo = log_choose(total, alt)

    lp_het = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_homalt = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)

    return lp_het, lp_homalt


def get_gt_metric2(events, mode, add_gt=True):

    params = {"DEL,TRA": [[0.5, 0.9], [0.2, 1 / 3.0]],
              "INS,DUP,INV,BND": [[0.4, 0.9], [0.6, 0.2]]}

    pp = {}
    for k, v in params.items():
        for kk in k.split(","):
            pp[kk] = v

    if add_gt:
        pass
    else:
        for e in events:
            e.GQ = '.'
            # e.SQ = '.'
            e.GT = './.'
        return events

    for e in events:

        a_params, b_params = pp[e.svtype]

        if e.svtype == "DEL" or e.svtype == "TRA":
            ref, support_reads = del_like(e)
        else:
            if mode == "pe":
                ref, support_reads = ins_like(e)
            else:
                ref, support_reads = ins_like_non_pe(e)

        if ref + support_reads <= 0:
            e.GQ = '.'
            # e.SQ = '.'
            e.GT = './.'
            continue

        gt_lplist = bayes_gt2(ref, support_reads, e.svtype == "DUP", a_params, b_params)
        best, second_best = sorted([(i, j) for i, j in enumerate(gt_lplist)], key=lambda x: x[1], reverse=True)[0:2]
        gt_idx = best[0]

        gt_sum = 0

        for gt in gt_lplist:
            try:
                gt_sum += 10**gt
            except OverflowError:
                gt_sum += 0
        if gt_sum > 0:
            gt_sum_log = math.log(gt_sum, 10)
            sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log))  # phred-scaled probability site is non-reference in this sample
            phred_gq = min(-10 * (second_best[1] - best[1]), 200)
            e.GQ = int(phred_gq)
            # e.SQ = sample_qual
            if gt_idx == 0:
                e.GT = '0/1'
            elif gt_idx == 1:
                e.GT = '1/1'
        else:
            e.GQ = '.'
            # e.SQ = '.'
            e.GT = './.'

    return events


def compressability(events):
    cdef EventResult_t e
    for e in events:
        c1 = []
        if e.contig:
            cont = e.contig.upper().encode("ascii")
            b = bytes(cont)
            c1.append(len(zlib.compress(b)) / len(b))
        if e.contig2:
            cont = e.contig2.upper().encode("ascii")
            b = bytes(cont)
            c1.append(len(zlib.compress(b)) / len(b))
        if c1:
            e.compress = round((sum(c1) / len(c1)) * 100, 2)
        else:
            e.compress = 0
    return events


def apply_model(df, mode, contigs, diploid, paired, thresholds):

    pth = os.path.dirname(os.path.abspath(__file__))
    pth = f"{pth}/dysgu_model.1.pkl.gz"
    models = pickle.load(gzip.open(pth, "rb"))

    if diploid == 'False' and contigs == 'False':
        raise NotImplemented("Choose either diploid == False or contigs == False, not both")

    # i.e. key is "nanopore_classifier_no_contigs"
    col_key = f"{mode}_cols{'_no_contigs' if contigs == 'False' else ''}{'_nodip' if diploid == 'False' else ''}"
    model_key = f"{mode}_classifier{'_no_contigs' if contigs == 'False' else ''}{'_nodip' if diploid == 'False' else ''}"

    cols = models[col_key]
    clf = models[model_key]
    logging.info(f"Model: {mode}, diploid: {diploid}, contig features: {contigs}. N features: {len(cols)}")

    c = dict(zip(
        ['NMS',    'SQC', 'CIPOS95',  'CIEND95',  'GC', 'REP', 'REPSC',  'SU', 'WR',       'SR',   'SC', 'NEXP',        'RPOLY',          'STRIDE', 'SVTYPE', 'SVLEN', 'NMP',   'NMB',    'MAPQP',   'MAPQS',    'NP', 'MAS',       'BE',         'COV',            'MCOV', 'NEIGH', 'NEIGH10',   'RB',        'PS',   'MS',    'NG',     'NSA',  'NXA',  'NMU',              'NDC',          'RMS',         'RED',      'BCC',            'STL',          'BND', 'SCW', 'RAS', 'FAS', 'OL',            'FCC', 'CMP',     'NG',        'RR',      'JIT',    'SBT',            'SQR'],
        ['NMsupp', 'sqc', 'cipos95A', 'cipos95B', 'gc', 'rep', 'rep_sc', 'su', 'spanning', 'supp', 'sc', 'n_expansion', 'ref_poly_bases', 'stride', 'svtype', 'svlen', 'NMpri', 'NMbase', 'MAPQpri', 'MAPQsupp', 'NP', 'maxASsupp', 'block_edge', 'raw_reads_10kb', 'mcov', 'neigh', 'neigh10kb', 'ref_bases', 'plus', 'minus', 'n_gaps', 'n_sa', 'n_xa', 'n_unmapped_mates', 'double_clips', 'remap_score', 'remap_ed', 'bad_clip_count', 'n_small_tlen', 'bnd', 'scw', 'ras', 'fas', "query_overlap", 'fcc', 'compress', "n_in_grp", 'ref_rep', 'jitter', 'strand_binom_t', 'clip_qual_ratio']
                 ))

    X = df[[c[i] for i in cols]]
    X.columns = cols
    keys = {"DEL": 1, "INS": 2, "DUP": 3, "INV": 4, "TRA": 2, "INV:DUP": 2, "BND": 4}

    X["SVTYPE"] = [keys[i] for i in X["SVTYPE"]]
    X["SVLEN"] = [i if i == i and i is not None else -1 for i in X["SVLEN"]]

    for c in models["cats"]:  # categorical data
        if c in X:
            X[c] = [i if i == i and i is not None else 0 for i in X[c]]
            X[c] = X[c].astype("category")

    pred = np.round(models[model_key].predict_proba(X)[:, 1], 3)
    df = df.assign(prob=pred)
    df = df.assign(filter=["PASS" if ((svt in thresholds and i >= thresholds[svt]) or (svt not in thresholds and i >= 0.5)) else "lowProb" for svt, i in zip(df["svtype"], pred)])

    return df


def bayes_multiple_observations(obs, prior=0.5):
    # goal is to calculate P(A | B1 & B2) where B1 and B2 are the two calls of a SV
    # calculate the event probability given the event was called twice
    # prior probability value of 0.5 as the classifier is assumed to be unbiased
    p_A_given_B = 0
    if len(obs) == 0:
        raise ValueError("post_call_metrics.pyx: len(obs) == 0")

    for p in obs:
        p_B = (p * prior) + ((1 - p) * (1-prior))
        p_A_given_B = (p * prior) / p_B
        prior = p_A_given_B
    return round(p_A_given_B, 3)


def update_prob_at_sites(df, events, thresholds, parse_probs, default_prob):
    # update probs for --sites using --sites-p or PROB/MeanPROB
    # allow only 1 call per matching variant in --sites
    new_probs = []
    best_idx = {}  # best df idex for each of matching --sites (use prob from this run to determine, not combined prob)
    for i, p, e in zip(df.index, df.prob, events):

        if e.site_info:
            site_id = e.site_info.id
            if site_id not in best_idx:
                best_idx[site_id] = (i, p)
            elif best_idx[site_id][1] < p:
                best_idx[site_id] = (i, p)

            if e.site_info.svtype == e.svtype:
                other_p = e.site_info.prob if parse_probs else default_prob
                obs = (p, other_p)
                p2 = bayes_multiple_observations(obs)
                new_probs.append(p2)
            else:
                new_probs.append(p)
        else:
            new_probs.append(p)

    df.prob = new_probs
    PASS = ['lowProb'] * len(df)
    for idx, (svt, i) in enumerate(zip(df["svtype"], new_probs)):
        if (svt in thresholds and i >= thresholds[svt]) or (svt not in thresholds and i >= 0.5):
            PASS.append('PASS')

    # df["PASS"] = ["PASS" if ((svt in thresholds and i >= thresholds[svt]) or (svt not in thresholds and i >= 0.5)) else "lowProb" for svt, i in zip(df["svtype"], new_probs)]

    keep = [True] * len(df)
    for i, e in enumerate(events):
        if e.site_info:
            best = best_idx[e.site_info.id]
            if best[0] != i:
                keep[i] = False
    df = df[keep]

    return df
