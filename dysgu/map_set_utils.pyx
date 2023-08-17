#cython: language_level=3

import click
import numpy as np
cimport numpy as np
import cython
import time
import logging
from libcpp.vector cimport vector as cpp_vector
from libcpp.pair cimport pair as cpp_pair

from libc.stdlib cimport abs as c_abs
from libc.math cimport fabs as c_fabs

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK
from libc.stdint cimport uint32_t, uint16_t, int16_t, int32_t

import math

ctypedef cpp_pair[int, int] cpp_item
ctypedef cpp_pair[long, int] cpp_long_item


def echo(*args):
    click.echo(args, err=True)


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        logging.info('%r  %2.2f ms' % (method.__name__, (te - ts) * 1000))
        return result
    return timed


def merge_intervals(intervals, srt=True, pad=0, add_indexes=False):
    """
    Merge a list of intervals, the expected format is a 3-tuple e.g. (chromosome, start, end). If add_indexes is
    set to True, merge_intervals expects a 4-tuple with the last item corresponding to an index variable

    :param intervals: The list of intervals to merge
    :type intervals: iterable
    :param srt: Sort the intervals by chromosome and start position
    :type srt: bool
    :param pad: Add a padding to intervals before merging. E.g. pad=10 subtracts 10 from start and adds 10 to end of interval
    :type pad: int
    :param add_indexes: Add the indexes of merged intervals to the output
    :type add_indexes: bool
    :return: list of merged intervals
    :rtype: list

    >>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [['chr1', 1, 5], ['chr2', 3, 5]]

    >>> a = [("chr1", 1, 10, 0), ("chr1", 9, 11, 1), ("chr1", 20, 30, 2)]
    >>> merge_intervals(a, add_indexes=True)
    >>> [('chr1', 1, 11, [0, 1]), ['chr1', 20, 30, [2]]]
    """
    if srt:
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))  # by chrom, start, end (index)
    else:
        sorted_by_lower_bound = intervals

    if pad:
        if not add_indexes:
            sorted_by_lower_bound = [[c, 0 if i - pad < 0 else i - pad, j + pad] for c, i, j in sorted_by_lower_bound]
        else:
            sorted_by_lower_bound = [[c, 0 if i - pad < 0 else i - pad, j + pad, k] for c, i, j, k in sorted_by_lower_bound]

    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            if not add_indexes:
                merged.append(higher)
            else:
                merged.append(list(higher)[:3] + [[higher[3]]])
            continue
        elif higher[0] != merged[-1][0]:  # Dont merge intervals on different chroms
            if not add_indexes:
                merged.append(higher)
            else:
                merged.append(list(higher)[:3] + [[higher[3]]])
        else:
            lower = merged[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= lower[2]:
                if not add_indexes:
                    merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]))
                else:
                    merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]), lower[3] + [higher[3]])
            else:
                if not add_indexes:
                    merged.append(higher)
                else:
                    merged.append(list(higher)[:3] + [[higher[3]]])
    return merged


cdef class Py_BasicIntervalTree:
    def __cinit__(self):
        self.thisptr = new BasicIntervalTree()
    def __dealloc__(self):
        del self.thisptr
    cpdef void add(self, int start, int end, int index):
        self.thisptr.add(start, end, index)
    cpdef bint searchInterval(self, int pos, int pos2):
        return self.thisptr.searchInterval(pos, pos2)
    cpdef overlappingInterval(self, int pos, int pos2):
        cdef Interval* res = self.thisptr.overlappingInterval(pos, pos2)
        if res[0] is None:  # [0] dereferences pointer
            return None
        else:
            return res[0].low, res[0].high
    cpdef void index(self):
        self.thisptr.index()
    cpdef allOverlappingIntervals(self, int start, int end):
        cdef cpp_vector[int] res
        self.thisptr.allOverlappingIntervals(start, end, res)
        return list(res)
    cpdef int countOverlappingIntervals(self, int pos, int pos2):
        return self.thisptr.countOverlappingIntervals(pos, pos2)


cdef class Py_DiGraph:
    """DiGraph, weighted"""
    def __cinit__(self):
        self.thisptr = new DiGraph()
    def __dealloc__(self):
        del self.thisptr
    cdef int addNode(self):
        return self.thisptr.addNode()
    cdef int hasEdge(self, int u, int v):
        return self.thisptr.hasEdge(u, v)
    cdef void addEdge(self, int u, int v, int w):
        self.thisptr.addEdge(u, v, w)
    cdef void updateEdge(self, int u, int v, int w):
        self.thisptr.addEdge(u, v, w)
    cdef int numberOfNodes(self) nogil:
        return self.thisptr.numberOfNodes()
    cdef void forInEdgesOf(self, int u, cpp_vector[cpp_pair[int, int]]& inEdges) nogil:
        self.thisptr.forInEdgesOf(u, inEdges)
    cdef void neighbors(self, int u, cpp_vector[int]& neigh) nogil:
        self.thisptr.neighbors(u, neigh)
    cdef float node_path_quality(self, int u, int v, int w)  nogil:
        return self.thisptr.node_path_quality(u, v, w)


cdef class Py_SimpleGraph:
    """Graph, weighted"""
    def __cinit__(self):
        self.thisptr = new SimpleGraph()
    def __dealloc__(self):
        del self.thisptr
    cpdef int addNode(self):
        return self.thisptr.addNode()
    cpdef int hasEdge(self, int u, int v):
        return self.thisptr.hasEdge(u, v)
    cpdef void addEdge(self, int u, int v, int w):
        self.thisptr.addEdge(u, v, w)
    cpdef int edgeCount(self):
        return self.thisptr.edgeCount()
    cpdef int weight(self, int u, int v):
        return self.thisptr.weight(u, v)
    cdef void neighbors(self, int u, cpp_vector[int]& neigh):
        self.thisptr.neighbors(u, neigh)
    cpdef void removeNode(self, int u):
        self.thisptr.removeNode(u)
    cdef void connectedComponents(self, char* pth, bint low_mem, cpp_vector[int]& components):
        self.thisptr.connectedComponents(pth, low_mem, components)
    cpdef int showSize(self):
        return self.thisptr.showSize()


cdef class Py_Int2IntMap:
    """Fast integer to integer unordered map using robin_hood flat map"""
    def __cinit__(self):
        self.thisptr = new Int2IntMap()
    def __dealloc__(self):
        del self.thisptr
    cdef void insert(self, int key, int value) nogil:
        self.thisptr.insert(key, value)
    cdef void erase(self, int key) nogil:
        self.thisptr.erase(key)
    cdef int has_key(self, int key) nogil:
        return self.thisptr.has_key(key)
    cdef int get(self, int key) nogil:
        return self.thisptr.get(key)
    cdef get_val_result get_value(self, int key) nogil:
        return self.thisptr.get_value(key)
    cdef int size(self) nogil:
        return self.thisptr.size()


cdef class Py_IntSet:
    """Fast set using robin_hood unordered set"""
    def __cinit__(self):
        self.thisptr = new IntSet()
    def __dealloc__(self):
        del self.thisptr
    cdef void insert(self, int key) nogil:
        self.thisptr.insert(key)
    cdef void erase(self, int key) nogil:
        self.thisptr.erase(key)
    cdef int has_key(self, int key) nogil:
        return self.thisptr.has_key(key)
    cdef int size(self) nogil:
        return self.thisptr.size()


cdef int cigar_exists(r):
    if r.cigartuples:
        return 1
    return 0


cdef tuple clip_sizes(r):
    c = r.cigartuples
    if not c:
        return 0, 0

    cdef int left = 0
    cdef int right = 0

    if c[0][0] == 4:
        left = c[0][1]
    if c[-1][0] == 4:
        right = c[-1][1]
    return left, right


cdef tuple clip_sizes_hard(r):
    c = r.cigartuples
    if not c:
        return 0, 0

    cdef int left = 0
    cdef int right = 0
    c1 = c[0][0]
    if c1 == 4 or c1 == 5:
        left = c[0][1]
    c1 = c[-1][0]
    if c1 == 4 or c1 == 5:
        right = c[-1][1]
    return left, right


cdef int cigar_clip(r, int clip_length):

    c = r.cigartuples
    if not c:
        return 0
    if (c[0][0] == 4 and c[0][1] >= clip_length) or (c[-1][0] == 4 and c[-1][1] >= clip_length):
        return 1
    return 0


cpdef int is_overlapping(int x1, int x2, int y1, int y2) nogil:
    return int(max(x1, y1) <= min(x2, y2))


cdef float min_fractional_overlapping(int x1, int x2, int y1, int y2):
    cdef int temp_v
    if x1 == x2 or y1 == y2:
        return 0
    if x2 < x1:
        temp_v = x2
        x2 = x1
        x1 = temp_v
    if y2 < y1:
        temp_v = y2
        y2 = y1
        y1 = temp_v
    cdef float overlap = float(max(0, (min(x2, y2) - max(x1, y1))))
    return min( overlap / float(c_abs(x2 - x1)),  overlap / float(c_abs(y2 - y1)) )

cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil:
    # Insertions have same x1/y1 position, use another measure
    if x1 == x2 or y1 == y2:
        return True

    cdef int temp_v
    if x2 < x1:
        temp_v = x2
        x2 = x1
        x1 = temp_v
    if y2 < y1:
        temp_v = y2
        y2 = y1
        y1 = temp_v
    cdef float overlap = float(max(0, (min(x2, y2) - max(x1, y1))))
    if overlap == 0:
        return False
    if (overlap / float(c_abs(x2 - x1))) > 0.1 and (overlap / float(c_abs(y2 - y1))) > 0.1:
        return True


cdef bint span_position_distance2(int x1, int x2, int y1, int y2):
    cdef int span1, span2, max_span
    cdef float span_distance, position_distance, center1, center2
    if x1 == x2:
        span1 = 1
        center1 = x1
    else:
        span1 = c_abs(x2 - x1)
        center1 = (x1 + x2) / 2
    if y1 == y2:
        span2 = 1
        center2 = y2
    else:
        span2 = c_abs(y2 - y1)
        center2 = (y1 + y2) / 2

    position_distance = c_fabs(center1 - center2)
    if position_distance > 2000:
        return 0
    max_span = max(span1, span2)
    span_distance = <float>c_abs(span1 - span2) / max_span
    if (position_distance / max_span) < 0.2 and span_distance < 0.3:
        return 1
    return 0


cdef bint span_position_distance(int x1, int x2, int y1, int y2, float norm, float thresh, ReadEnum_t read_enum,
                                 bint paired_end, int cigar_len1, int cigar_len2, bint trust_ins_len) nogil:
    # https://github.com/eldariont/svim/blob/master/src/svim/SVIM_clustering.py
    cdef int span1, span2, max_span
    cdef float span_distance, position_distance, center1, center2
    # if read_enum == BREAKEND:
    #     return 0
    cdef bint within_read_event = cigar_len1 > 0 and cigar_len2 > 0
    if within_read_event and c_abs(cigar_len1 - cigar_len2) > 500:
        return 0
    if read_enum == INSERTION and within_read_event and trust_ins_len:
        span1 = cigar_len1
        span2 = cigar_len2
        center1 = x1
        center2 = y1

    else:
        if x1 == x2:
            span1 = 1
            center1 = x1
        else:
            span1 = c_abs(x2 - x1)
            center1 = (x1 + x2) / 2
        if y1 == y2:
            span2 = 1
            center2 = y2
        else:
            span2 = c_abs(y2 - y1)
            center2 = (y1 + y2) / 2

    max_span = max(span1, span2)
    position_distance = c_fabs(center1 - center2)
    if position_distance > 2000:
        return 0

    span_distance = <float>c_abs(span1 - span2) / max_span

    if not paired_end or read_enum == SPLIT:
        # span position distance:
        if (position_distance / norm) + span_distance < thresh:
            return 1
        return 0

    if (position_distance / max_span) < thresh and span_distance < thresh:
        return 1
    return 0


cdef float position_distance(int x1, int x2, int y1, int y2) nogil:
    # https://github.com/eldariont/svim/blob/master/src/svim/SVIM_clustering.py
    cdef int span1, span2
    cdef float center1, center2
    if x1 == x2:
        span1 = 1
        center1 = x1
    else:
        span1 = c_abs(x2 - x1)
        center1 = (x1 + x2) / 2
    if y1 == y2:
        span2 = 1
        center2 = y2
    else:
        span2 = c_abs(y2 - y1)
        center2 = (y1 + y2) / 2
    return c_fabs(center1 - center2)


def to_dict(self):
    return {v: self.__getattribute__(v) for v in dir(self) if "__" not in v and v != "to_dict" and v != "from_dict"}


def from_dict(self, d):
    allowed = set(dir(self))
    for k, v in d.items():
        if k in allowed:
            self.__setattr__(k, v)
    return self


@cython.auto_pickle(True)
cdef class EventResult:
    """Data holder for classifying alignments into SV types"""
    def __init__(self):
        # set a few variables
        self.svlen_precise = 1
        self.rep = 0
        self.ref_rep = 0
        self.NMpri = 0
        self.NMbase = 0
        self.mcov = 0
        self.neigh = 0
        self.ref_bases = 0
        self.n_sa = 0
        self.n_gaps = 0
        self.compress = 0

    def __repr__(self):
        return str(to_dict(self))
        # return str(self.to_dict())

    # def __getstate__(self):  # for pickling
    #     return to_dict(self)
    #     # return self.to_dict()
    #
    # def __setstate__(self, d):
    #     for k, v in d.items():
    #         self.__setattr__(k, v)
