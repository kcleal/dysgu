#cython: language_level=3

import click
import numpy as np
cimport numpy as np
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


cdef class Py_CoverageTrack:

    def __init__(self, outpath, infile, max_coverage):
        self.current_chrom = -1
        self.outpath = outpath
        self.cov_array = np.array([], dtype="int32")
        self.infile = infile
        self.max_coverage = max_coverage

    def add(self, a):
        if a.flag & 1284 or a.cigartuples is None or a.mapq == 0:  # not primary, duplicate or unmapped?
            return

        if self.current_chrom != a.rname:
            self.write_track()
            self.set_cov_array(a.rname)
            self.current_chrom = a.rname

        cdef int32_t[:] arr = self.cov_array
        if int(a.pos  / 10) > len(arr) - 1:
            return

        cdef int index_start = a.pos
        cdef int index_bin = int(index_start / 10)
        cdef int opp, length

        for opp, length in a.cigartuples:
            if index_bin > len(arr) - 1:
                break
            if opp == 2:
                index_start += length
                index_bin = int(index_start / 10)
            elif opp == 0 or opp == 7 or opp == 8:
                arr[index_bin] += 1
                index_start += length
                index_bin = int(index_start / 10)
                if index_bin < len(arr):
                    arr[index_bin] -= 1

    def set_cov_array(self, int rname):
        chrom_length = self.infile.get_reference_length(self.infile.get_reference_name(rname))
        self.cov_array = np.resize(self.cov_array, int(chrom_length / 10) + 1)
        self.cov_array.fill(0)

    def write_track(self):
        cdef np.ndarray[int16_t, ndim=1] ca = self.cov_array.astype("int16")  # old style buffer for .tofile function
        cdef int current_cov = 0
        cdef int16_t v
        cdef int i
        if self.current_chrom != -1:
            out_path = "{}/{}.dysgu_chrom.bin".format(self.outpath, self.infile.get_reference_name(self.current_chrom))
            for i in range(len(ca)):
                v = ca[i]
                current_cov += v
                if current_cov < 32000:
                    ca[i] = current_cov
                elif current_cov < 0:
                    ca[i] = 0  # sanity check
                elif current_cov >= 32000:
                    ca[i] = 32000
            ca.tofile(out_path)


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
    cdef cpp_vector[cpp_pair[int, int]] forInEdgesOf(self, int u) nogil:
        return self.thisptr.forInEdgesOf(u)
    cdef cpp_vector[int] neighbors(self, int u) nogil:
        return self.thisptr.neighbors(u)
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
    cpdef cpp_vector[int] neighbors(self, int u):
        return self.thisptr.neighbors(u)
    cpdef void removeNode(self, int u):
        self.thisptr.removeNode(u)
    cpdef cpp_vector[int] connectedComponents(self):
        return self.thisptr.connectedComponents()
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


cdef bint span_position_distance2(int x1, int x2, int y1, int y2): # nogil:
    # https://github.com/eldariont/svim/blob/master/src/svim/SVIM_clustering.py
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

    position_distance = c_fabs(center1 - center2) # 1 #distance_normalizer
    if position_distance > 2000:
        return 0
    max_span = max(span1, span2)
    span_distance = <float>c_abs(span1 - span2) / max_span
    # echo("pd", position_distance, center1, center2)
    # echo((position_distance / max_span), span_distance, center1, center2 )
    if (position_distance / max_span) < 0.2 and span_distance < 0.3:
    # if position_distance < 100 and span_distance < 0.08:
        return 1
    return 0


cdef bint span_position_distance(int x1, int x2, int y1, int y2, float norm, float thresh, ReadEnum_t read_enum,
                                 bint paired_end) nogil:
    # https://github.com/eldariont/svim/blob/master/src/svim/SVIM_clustering.py
    cdef int span1, span2, max_span
    cdef float span_distance, position_distance, center1, center2
    if read_enum == BREAKEND:
        return 0
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

    if (position_distance / max_span) < thresh and span_distance < thresh:  # 0.2, 0.3
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


cdef class EventResult:
    """Data holder for classifying alignments into SV types"""
    def __cinit__(self):
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

    def to_dict(self):
        return {v: self.__getattribute__(v) for v in dir(self) if "__" not in v and v != "to_dict"}

    def from_dict(self, d):
        allowed = set(dir(self))
        for k, v in d.items():
            if k in allowed:
                self.__setattr__(k, v)
        return self

    def __repr__(self):
        return str(self.to_dict())
