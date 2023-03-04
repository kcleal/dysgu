#cython: language_level=3, boundscheck=True, c_string_type=unicode, c_string_encoding=utf8, infer_types=True
import numpy as np
cimport numpy as np
from dysgu.map_set_utils cimport unordered_map, EventResult
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement
from libc.math cimport fabs as c_fabs
from libcpp.vector cimport vector as cpp_vector
from dysgu.map_set_utils import echo
import array
import pandas as pd
pd.options.mode.chained_assignment = None
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

