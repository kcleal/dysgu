#cython: language_level=3

import click

from libcpp.vector cimport vector as cpp_vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.string cimport string as cpp_string

import cython

from libc.stdlib cimport abs as c_abs
from libc.math cimport fabs as c_fabs

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK
from libc.stdint cimport uint32_t


# ctypedef cpp_vector[int] int_vec_t
# ctypedef cpp_pair[int, int] get_val_result

# ctypedef Py_Int2IntVecMap[int, int_vec_t] node_dict_t
# ctypedef Py_IntVec2IntMap[int_vec_t, int] node_dict2_r_t


def echo(*args):
    click.echo(args, err=True)


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


# cdef class Py_StrSet:
#
#     def __cinit__(self):
#         self.thisptr = new StrSet()
#     def __dealloc__(self):
#         del self.thisptr
#     cdef void insert(self, cpp_string key) nogil:
#         self.thisptr.insert(key)
#     cdef void erase(self, cpp_string key) nogil:
#         self.thisptr.erase(key)
#     cdef int has_key(self, cpp_string key) nogil:
#         return self.thisptr.has_key(key)
#     cdef int size(self) nogil:
#         return self.thisptr.size()


# cdef class Py_PairScope:
#     """Cluster reads by their mate position"""
#     def __cinit__(self):
#         self.thisptr = new PairScope()
#     def __dealloc__(self):
#         del self.thisptr
#     cpdef void add_params(self, int m, int n):
#         self.thisptr.add_params(m, n)
#     cpdef cpp_vector[int] update(self, int node_name, int c_chrom, int c_pos, int chrom2, int pos2):
#         return self.thisptr.update(node_name, c_chrom, c_pos, chrom2, pos2)



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


cdef int cigar_clip(r, int clip_length):

    c = r.cigartuples
    if not c:
        return 0
    if (c[0][0] == 4 and c[0][1] >= clip_length) or (c[-1][0] == 4 and c[-1][1] >= clip_length):
        return 1
    return 0


cdef int is_overlapping(int x1, int x2, int y1, int y2) nogil:
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


cdef bint span_position_distance(int x1, int x2, int y1, int y2) nogil:
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
