#cython: language_level=3, c_string_type=unicode, c_string_encoding=utf8
import numpy as np
cimport numpy as np
from dysgu.map_set_utils cimport unordered_map, EventResult, cigar_clip, clip_sizes
from dysgu.map_set_utils import merge_intervals
from dysgu.io_funcs import intersecter #, iitree
from superintervals import IntervalSet
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement
from libc.math cimport fabs as c_fabs
from libc.stdint cimport uint32_t
from libcpp.vector cimport vector as cpp_vector
from dysgu.map_set_utils import echo
import array
from collections import defaultdict
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_cigar
import pandas as pd
pd.options.mode.chained_assignment = None
from sys import byteorder
import os
import resource
ctypedef EventResult EventResult_t


class BadClipCounter:

    def __init__(self, n_references, low_mem, temp_dir):
        self.low_mem = low_mem
        self.temp_dir = temp_dir
        self.n_refs = n_references
        if not self.low_mem:
            self.clip_pos_arr = [array.array("L", []) for i in range(n_references)]  # 32 bit unsigned int
        else:
            soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
            if soft + (n_references * 2) > soft:
                resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, soft + (n_references * 2)), hard))
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
        if chrom < 0 or len_a == 0:
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


cdef float soft_clip_qual_corr(reads):
    """Function to compute the correlation of quality values between the soft-clipped portions of reads. values of
    1 imply no correlation whereas closer to 0 implies that quality values are correlated"""
    cdef const unsigned char[:] quals
    cdef int idx, x, i, j
    cdef unordered_map[int, cpp_vector[int]] qq
    cdef int left_clip, right_clip
    cdef AlignedSegment r
    for r in reads:
        quals = r.query_qualities
        left_clip = 0
        right_clip = 0
        clip_sizes(r, &left_clip, &right_clip)

        if quals is None:
            continue
        if left_clip:
            idx = 0
            for x in range(r.pos - left_clip, r.pos):
                qq[x].push_back(quals[idx])
                idx += 1
        if right_clip:
            end_pos = r.query_alignment_end + right_clip
            for idx in range(-1, - right_clip, -1):
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
    cdef int size;
    with nogil:
        qq_iter = qq.begin()
        while qq_iter != qq.end():
            second = dereference(qq_iter).second
            size = second.size()
            if size > 1:
                sum_second = 0
                for i in range(size):
                    sum_second += second[i]
                mean_second = sum_second / float(size)
                for i in range(size):
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


cdef int CMATCH = 0
cdef int CINS = 1
cdef int CDEL = 2
cdef int CSOFT_CLIP = 4
cdef int CHARD_CLIP = 5
cdef int CEQUAL = 7
cdef int CDIFF = 8

cdef int BAM_CIGAR_SHIFT = 4
cdef int BAM_CIGAR_MASK  = 15 #0xf


# https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
cdef int[32] tab32 = [
    0,  9,  1, 10, 13, 21,  2, 29,
    11, 14, 16, 18, 22, 25,  3, 30,
     8, 12, 20, 28, 15, 17, 24,  7,
    19, 27, 23,  6, 26,  5,  4, 31]

cdef int log2_32(uint32_t value):
    value |= value >> 1
    value |= value >> 2
    value |= value >> 4
    value |= value >> 8
    value |= value >> 16
    return tab32[<uint32_t>(value*0x07C4ACDD) >> 27]


cdef void window_rate(WindowRate *result, uint32_t cigar_l, uint32_t *cigar_p, int index, int window_size, bint reverse):
    cdef int n = 0
    cdef int matches = 0
    cdef int covered = 0
    cdef int i = index
    cdef int stop, opp, l
    cdef uint32_t cigar_value
    if reverse:
        stop = 0
    else:
        stop = cigar_l - 1
    while i != stop:
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & BAM_CIGAR_MASK
        l = <int> cigar_value >> BAM_CIGAR_SHIFT
        if opp == CMATCH or opp == CEQUAL:
            matches += l
            covered += l
        elif opp == CDEL:
            covered += l
            n += 1
        elif opp == CINS or opp == CSOFT_CLIP or opp == CHARD_CLIP or opp == CDIFF:
            n += 1
        if reverse:
            i -= 1
        else:
            i += 1
        if covered >= window_size:
            # if final match block is really long, the previous cigar opp can be missed if this is not done:
            matches = min(window_size, matches)
            break
    result.rate = 0 if not matches else n / (matches + n)
    result.index = i


cpdef filter_poorly_aligned_ends(spanning_alignments, float divergence=0.02):
    # This trys to check if the ends of reads are poorly aligned by looking at the rate of non-matching cigar events
    cdef int i, index_begin, index_end, start_i, end_i, cigar_index
    cdef float threshold, w_rate, mean, std
    cdef uint32_t cigar_value
    cdef uint32_t cigar_l
    cdef uint32_t *cigar_p
    cdef WindowRate window_r
    cdef AlignedSegment alignment
    rates = []
    for item in spanning_alignments:
        alignment = item.align
        if alignment.flag & 2048:
            continue
        cigar_l = alignment._delegate.core.n_cigar
        cigar_p = bam_get_cigar(alignment._delegate)
        window_rate(&window_r, cigar_l, cigar_p, 0, 1_000_000_000, <bint>False)
        rates.append(window_r.rate)
    if rates:
        mean = np.mean(rates)
        std = np.std(rates)
    else:
        mean = 0
        std = 0

    multiplier = 3
    threshold = min(divergence, max(0.008, multiplier*mean, mean + (std*3)))
    spanning = []
    for item in spanning_alignments:
        alignment = item.align
        cigar_l = alignment._delegate.core.n_cigar
        cigar_p = bam_get_cigar(alignment._delegate)
        cigar_index = item.cigar_index
        index_begin = 0
        for i in range(cigar_l):
            window_rate(&window_r, cigar_l, cigar_p, i, 2000, <bint>False)
            w_rate = window_r.rate
            start_i = window_r.index
            if w_rate < threshold:
                break
            index_begin = start_i
        index_end = cigar_l
        for i in range(cigar_l-1, index_begin -1, -1):
            window_rate(&window_r, cigar_l, cigar_p, i, 2000, <bint>True)
            w_rate = window_r.rate
            end_i = window_r.index
            if w_rate < threshold:
                break
            index_end = end_i
        if index_begin <= cigar_index <= index_end:
            spanning.append(item)
    return spanning, 1 - (len(spanning) / len(spanning_alignments))


cpdef gap_size_upper_bound(AlignedSegment alignment, int cigarindex, int pos_input, int end_input, int length_extend=15, float divergence=0.02, float len_t=3):
    # expand indel using cigar, merges nearby gaps into the SV event
    cdef int pos, end, l, opp, extent_left, extent_right, candidate_type, candidate_len, len_input, i, dist, dist_thresh, middle, last_seen_size
    cdef uint32_t cigar_value
    cdef uint32_t cigar_l = alignment._delegate.core.n_cigar
    cdef uint32_t *cigar_p = bam_get_cigar(alignment._delegate)

    pos = pos_input
    end = end_input
    extent_left = pos
    extent_right = end
    cigar_value = cigar_p[cigarindex]
    candidate_type = cigar_value & BAM_CIGAR_MASK
    candidate_len = cigar_value >> BAM_CIGAR_SHIFT
    len_input = candidate_len
    dist_thresh = min(<int> (len_input * len_t), <int> (<float> log2_32(1 + candidate_len) / divergence))
    i = cigarindex + 1
    dist = 0
    last_seen_size = candidate_len

    while i < cigar_l:  # extent rightwards
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & BAM_CIGAR_MASK
        l = <int> cigar_value >> BAM_CIGAR_SHIFT
        if opp == CSOFT_CLIP or opp == CHARD_CLIP:
            break
        if opp == CMATCH or opp == CEQUAL or opp == CDIFF:
            dist += l
            end += l
            if dist > dist_thresh:
                break
        elif l > length_extend and (opp == CDEL or opp == CINS):
            if candidate_type == CDEL:
                if opp == CDEL:
                    candidate_len += l
                    end += l
                    extent_right = end
                else:
                    if l >= candidate_len: break
                    candidate_len -= l
            else:
                if opp == CDEL:
                    if l >= candidate_len: break
                    candidate_len -= l
                    end -= l
                else:
                    candidate_len += l
                    extent_right = end
            dist_thresh = min(<int> (l * len_t), <int> (<float> log2_32(1 + candidate_len) / divergence))
            dist = 0
        i += 1

    i = cigarindex - 1
    dist = 0
    while i > -1:  # extend leftwards
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & BAM_CIGAR_MASK
        l = <int> cigar_value >> BAM_CIGAR_SHIFT
        if opp == CSOFT_CLIP or opp == CHARD_CLIP:
            break
        if opp == CMATCH or opp == CEQUAL or opp == CDIFF:
            dist += l
            pos -= l
            if dist > dist_thresh:
                break
        elif l > length_extend and (opp == CDEL or opp == CINS):
            if candidate_type == CDEL:
                if opp == CDEL:
                    candidate_len += l
                    pos -= l
                    extent_left = pos
                else:
                    if l >= candidate_len: break
                    candidate_len -= l
            else:
                if opp == CDEL:
                    if l >= candidate_len: break
                    candidate_len -= l
                    pos -= l
                else:
                    candidate_len += l
                    extent_left = pos
            dist_thresh = min(<int> (l * len_t), <int>( <float>log2_32(1 + candidate_len) / divergence) )
            dist = 0
        i -= 1

    if extent_right > extent_left and (extent_left != pos_input or extent_right != end_input):  # has been extended
        # middle = extent_left + <int>((extent_right - extent_left) / 2)
        # if candidate_type == CDEL:
        #     pos = middle - int(candidate_len * 0.5)
        #     end = middle + int(candidate_len * 0.5)
        # else:
        #     pos = pos_input
        #     end = end_input
        # return candidate_len, pos, end
        return candidate_len, extent_left, extent_right

    return len_input, pos_input, end_input


def sample_level_density(potential, regions, max_dist=50):
    tmp_list = defaultdict(list)
    cdef EventResult_t ei
    cdef int start, stop, idx

    for idx in range(len(potential)):
        ei = potential[idx]
        # Only find density for non-region calls, otherwise too dense to be meaningful
        if not intersecter(regions, ei.chrA, ei.posA, ei.posA + 1):
            tmp_list[ei.chrA].append((ei.posA - max_dist, ei.posA + max_dist, idx))
        if not intersecter(regions, ei.chrB, ei.posB, ei.posB + 1):
            tmp_list[ei.chrB].append((ei.posB - max_dist, ei.posB + max_dist, idx))
    # nc2 = {k: iitree(v, add_value=True) for k, v in tmp_list.items()}

    si_sets = {}

    for k, v in tmp_list.items():
        iset = IntervalSet(with_data=True)
        for start, stop, idx in v:
            iset.add_int_value(start, stop, idx)
        iset.index()
        si_sets[k] = iset

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
            # vv = nc2[ei.chrA].countOverlappingIntervals(ei.posA, ei.posA + 1)
            vv = si_sets[ei.chrA].count_overlaps(ei.posA, ei.posA + 1) #countOverlappingIntervals(ei.posA, ei.posA + 1)
            neighbors += vv - expected
            count += 1
        if not intersecter(regions, ei.chrB, ei.posB, ei.posB + 1):
            # vv = nc2[ei.chrB].countOverlappingIntervals(ei.posB, ei.posB + 1)
            vv = si_sets[ei.chrB].count_overlaps(ei.posB, ei.posB + 1)
            neighbors += vv - expected
            count += 1
        neighbors_10kb = 0.
        count_10kb = 0
        large_itv = merge_intervals(((ei.chrA, ei.posA, ei.posA + 1), (ei.chrB, ei.posB, ei.posB + 1)), pad=10000)
        for c, s, e in large_itv:
            if not intersecter(regions, c, s, e):
                # vv = nc2[c].countOverlappingIntervals(s, e)
                vv = si_sets[c].count_overlaps(s, e)
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


cpdef dict search_ssr_kc(ori):
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
