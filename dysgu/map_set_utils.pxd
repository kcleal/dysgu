#cython: language_level=3

from libcpp.vector cimport vector as cpp_vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.string cimport string as cpp_string

import cython

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK


from libc.stdint cimport uint32_t, uint8_t, uint64_t

ctypedef cpp_vector[int] int_vec_t
ctypedef cpp_pair[int, int] get_val_result

# ctypedef Py_Int2IntVecMap[int, int_vec_t] node_dict_t
# ctypedef Py_IntVec2IntMap[int_vec_t, int] node_dict2_r_t

from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators


cdef extern from "xxhash64.h" namespace "XXHash64" nogil:
    #static uint64_t hash(const void* input, uint64_t length, uint64_t seed)
    cdef uint64_t hash(void* input, uint64_t length, uint64_t seed) nogil


cdef extern from "robin_hood.h" namespace "robin_hood" nogil:
    cdef cppclass unordered_set[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void insert(T&)
        void erase(T&)
        int size()
        iterator find(const T&)
        iterator begin()
        iterator end()
        void clear()
        bint empty()


# cdef extern from "robin_set.h" namespace "tsl" nogil:
#     cdef cppclass robin_set[T]:
#         cppclass iterator:
#             T operator*()
#             iterator operator++()
#             bint operator==(iterator)
#             bint operator!=(iterator)
#         vector()
#         void insert(T&)
#         void erase(T&)
#         int size()
#         iterator find(const T&, size_t precalculated_hash)
#         iterator begin()
#         iterator end()
#         void clear()
#         bint empty()


cdef extern from "wrap_map_set2.h" nogil:
    cdef cppclass DiGraph:
        DiGraph() nogil

        int addNode()
        int hasEdge(int, int)
        void addEdge(int, int, int)
        int weight(int, int)
        void updateEdge(int, int, int)
        int numberOfNodes() nogil
        cpp_vector[cpp_pair[int, int]] forInEdgesOf(int) nogil
        cpp_vector[int] neighbors(int) nogil
        float node_path_quality(int, int, int) nogil


cdef class Py_DiGraph:
    """DiGraph, no weight"""
    cdef DiGraph *thisptr

    cdef int addNode(self)
    cdef int hasEdge(self, int u, int v)
    cdef void addEdge(self, int u, int v, int w)
    cdef void updateEdge(self, int u, int v, int w)
    cdef int numberOfNodes(self) nogil
    cdef cpp_vector[cpp_pair[int, int]] forInEdgesOf(self, int u) nogil
    cdef cpp_vector[int] neighbors(self, int u) nogil
    cdef float node_path_quality(self, int u, int v, int w) nogil

cdef extern from "wrap_map_set2.h":
    cdef cppclass SimpleGraph:
        SimpleGraph()

        int addNode()
        int hasEdge(int, int)
        void addEdge(int, int, int)
        int edgeCount()
        int weight(int, int)
        cpp_vector[int] neighbors(int)
        void removeNode(int)
        cpp_vector[int] connectedComponents()
        int showSize()



cdef class Py_SimpleGraph:
    """Graph"""
    cdef SimpleGraph *thisptr

    cpdef int addNode(self)
    cpdef int hasEdge(self, int u, int v)
    cpdef void addEdge(self, int u, int v, int w)
    cpdef int edgeCount(self)
    cpdef int weight(self, int u, int v)
    cpdef cpp_vector[int] neighbors(self, int u)
    cpdef void removeNode(self, int u)
    cpdef cpp_vector[int] connectedComponents(self)
    cpdef int showSize(self)


cdef extern from "wrap_map_set2.h" nogil:
    cdef cppclass Int2IntMap:
        Int2IntMap() nogil
        void insert(int, int) nogil
        void erase(int) nogil
        int has_key(int) nogil
        int get(int) nogil
        get_val_result get_value(int key) nogil
        int size() nogil


cdef class Py_Int2IntMap:
    """Fast integer to integer unordered map using tsl::robin-map"""
    cdef Int2IntMap *thisptr

    cdef void insert(self, int key, int value) nogil
    cdef void erase(self, int key) nogil
    cdef int has_key(self, int key) nogil
    cdef int get(self, int key) nogil
    cdef get_val_result get_value(self, int key) nogil
    cdef int size(self) nogil


cdef extern from "wrap_map_set2.h" nogil:
    cdef cppclass IntSet:
        IntSet() nogil
        void insert(int) nogil
        void erase(int) nogil
        int has_key(int) nogil
        int get(int) nogil
        int size() nogil


cdef class Py_IntSet:
    """Fast 32 bit int set using tsl::robin-set"""
    cdef IntSet *thisptr

    cdef void insert(self, int key) nogil
    cdef void erase(self, int key) nogil
    cdef int has_key(self, int key) nogil
    cdef int size(self) nogil


# cdef extern from "wrap_map_set2.h":
#     cdef cppclass StrSet:
#         StrSet()
#         void insert(cpp_string) nogil
#         void erase(cpp_string) nogil
#         int has_key(cpp_string) nogil
#         int size() nogil
#
#
# cdef class Py_StrSet:
#     """Fast std::string set using tsl::robin-set"""
#     cdef StrSet *thisptr
#     cdef void insert(self, cpp_string key) nogil
#     cdef void erase(self, cpp_string key) nogil
#     cdef int has_key(self, cpp_string key) nogil
#     cdef int size(self) nogil


# cdef extern from "wrap_map_set2.h" nogil:
#     cdef cppclass PairScope:
#         PairScope() nogil
#         void add_params(int, int)
#         cpp_vector[int] update(int, int, int, int, int)
#
# cdef class Py_PairScope:
#     cpdef PairScope * thisptr
#     cpdef void add_params(self, int, int)
#     cpdef cpp_vector[int] update(self, int, int, int, int, int)


cdef int cigar_exists(r)


cdef tuple clip_sizes(r)


cdef tuple clip_sizes_hard(r)


cdef int cigar_clip(r, int clip_length)


cdef int is_overlapping(int x1, int x2, int y1, int y2) nogil


cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil


cdef bint span_position_distance(int x1, int x2, int y1, int y2) nogil


cdef float position_distance(int x1, int x2, int y1, int y2) nogil


# cdef list search_ssr_kc(str seq)

cdef void sliding_window_minimum(int k, int m, str s, unordered_set[long]& found)