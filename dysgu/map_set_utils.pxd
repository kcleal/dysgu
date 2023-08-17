#cython: language_level=3

from libcpp.vector cimport vector as cpp_vector
from libcpp.pair cimport pair as cpp_pair
from libcpp.utility cimport pair

import numpy as np
cimport numpy as np

from libc.stdint cimport uint64_t, int32_t, int8_t

ctypedef cpp_vector[int] int_vec_t
ctypedef cpp_pair[int, int] get_val_result


ctypedef enum ReadEnum_t:
    DISCORDANT = 0
    SPLIT = 1
    DELETION = 2
    INSERTION = 3
    BREAKEND = 4


cdef extern from "xxhash64.h" namespace "XXHash64" nogil:
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


cdef extern from "robin_hood.h" namespace "robin_hood" nogil:
    cdef cppclass unordered_map[T, U, HASH=*]:
        ctypedef T key_type
        ctypedef U mapped_type
        ctypedef pair[const T, U] value_type
        cppclass iterator:
            pair[T, U]& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)
        unordered_map() except +
        unordered_map(unordered_map&) except +
        U& operator[](T&)
        pair[iterator, bint] insert(pair[T, U])
        iterator find(const T&)
        iterator begin()
        iterator end()
        int erase(T&)
        int size()
        void clear()
        bint empty()


cdef extern from "find_reads.hpp" nogil:
    cdef cppclass CoverageTrack:
        CoverageTrack() nogil

        void add(int, int)
        int get_cov(int)
        bint cov_val_good(int, int, int)
        void set_cov_array(int)
        void write_track(char*)
        void set_max_cov(int)


cdef extern from "graph_objects.hpp":
    ctypedef struct Interval:
        int low
        int high


cdef extern from "graph_objects.hpp":
    cdef cppclass BasicIntervalTree:
        BasicIntervalTree()
        void add(int, int, int)
        bint searchInterval(int, int)
        Interval* overlappingInterval(int, int)
        void index()
        void allOverlappingIntervals(int, int, cpp_vector[int])
        int countOverlappingIntervals(int, int)


cdef class Py_BasicIntervalTree:
    cdef BasicIntervalTree *thisptr

    cpdef void add(self, int start, int end, int index)
    cpdef bint searchInterval(self, int pos, int pos2)
    cpdef overlappingInterval(self, int pos, int pos2)
    cpdef void index(self)
    cpdef allOverlappingIntervals(self, int pos, int pos2)
    cpdef int countOverlappingIntervals(self, int pos, int pos2)


cdef extern from "graph_objects.hpp" nogil:
    cdef cppclass DiGraph:
        DiGraph() nogil

        int addNode()
        int hasEdge(int, int)
        void addEdge(int, int, int)
        int weight(int, int)
        void updateEdge(int, int, int)
        int numberOfNodes() nogil
        void forInEdgesOf(int, cpp_vector[cpp_pair[int, int]]&) nogil
        void neighbors(int, cpp_vector[int]&) nogil
        float node_path_quality(int, int, int) nogil


cdef class Py_DiGraph:
    """DiGraph, no weight"""
    cdef DiGraph *thisptr

    cdef int addNode(self)
    cdef int hasEdge(self, int u, int v)
    cdef void addEdge(self, int u, int v, int w)
    cdef void updateEdge(self, int u, int v, int w)
    cdef int numberOfNodes(self) nogil
    cdef void forInEdgesOf(self, int u, cpp_vector[cpp_pair[int, int]]& inEdges) nogil
    cdef void neighbors(self, int u, cpp_vector[int]& neigh) nogil
    cdef float node_path_quality(self, int u, int v, int w) nogil


cdef extern from "graph_objects.hpp" nogil:
    cdef cppclass MinimizerTable:
        MinimizerTable() nogil

        int size()
        void insert(long key, long value1)
        void erase(long key)
        void erase_lower(long key, long value)
        int has_key(long key)
        int has_lower_key(long key2)
        long get_lower()
        unordered_set[long].iterator get_iterator()
        unordered_set[long].iterator get_iterator_begin()
        unordered_set[long].iterator get_iterator_end()


cdef extern from "graph_objects.hpp":
    cdef cppclass SimpleGraph:
        SimpleGraph()

        int addNode()
        int hasEdge(int, int)
        void addEdge(int, int, int)
        int edgeCount()
        int weight(int, int)
        void neighbors(int, cpp_vector[int]&)
        void removeNode(int)
        void connectedComponents(char*, bint, cpp_vector[int]&)
        int showSize()


cdef class Py_SimpleGraph:
    """Graph"""
    cdef SimpleGraph *thisptr

    cpdef int addNode(self)
    cpdef int hasEdge(self, int u, int v)
    cpdef void addEdge(self, int u, int v, int w)
    cpdef int edgeCount(self)
    cpdef int weight(self, int u, int v)
    cdef void neighbors(self, int u, cpp_vector[int]& neigh)
    cpdef void removeNode(self, int u)
    cdef void connectedComponents(self, char* pth, bint low_mem, cpp_vector[int]& neigh)
    cpdef int showSize(self)


cdef extern from "graph_objects.hpp" nogil:
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


cdef extern from "graph_objects.hpp" nogil:
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


cdef extern from "<map>" namespace "std" nogil:
    cdef cppclass multimap[T, U, COMPARE=*, ALLOCATOR=*]:
        ctypedef T key_type
        ctypedef U mapped_type
        ctypedef pair[const T, U] value_type
        ctypedef COMPARE key_compare
        ctypedef ALLOCATOR allocator_type
        cppclass iterator:
            pair[T, U]& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)
        cppclass reverse_iterator:
            pair[T, U]& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(reverse_iterator)
            bint operator!=(reverse_iterator)
        cppclass const_iterator(iterator):
            pass
        cppclass const_reverse_iterator(reverse_iterator):
            pass
        multimap() except +
        multimap(multimap&) except +
        U& operator[](T&)
        bint operator==(multimap&, multimap&)
        bint operator!=(multimap&, multimap&)
        bint operator<(multimap&, multimap&)
        bint operator>(multimap&, multimap&)
        bint operator<=(multimap&, multimap&)
        bint operator>=(multimap&, multimap&)
        U& at(const T&) except +
        const U& const_at "at"(const T&) except +
        iterator begin()
        const_iterator const_begin "begin" ()
        void clear()
        size_t count(const T&)
        bint empty()
        iterator end()
        const_iterator const_end "end" ()
        pair[iterator, iterator] equal_range(const T&)
        #pair[const_iterator, const_iterator] equal_range(key_type&)
        void erase(iterator)
        void erase(iterator, iterator)
        size_t erase(const T&)
        iterator find(const T&)
        const_iterator const_find "find" (const T&)
        pair[iterator, bint] insert(pair[T, U]) except + # XXX pair[T,U]&
        iterator insert(iterator, pair[T, U]) except + # XXX pair[T,U]&
        #void insert(input_iterator, input_iterator)
        #key_compare key_comp()
        iterator lower_bound(const T&)
        const_iterator const_lower_bound "lower_bound"(const T&)
        size_t max_size()
        reverse_iterator rbegin()
        const_reverse_iterator const_rbegin "rbegin"()
        reverse_iterator rend()
        const_reverse_iterator const_rend "rend"()
        size_t size()
        void swap(multimap&)
        iterator upper_bound(const T&)
        const_iterator const_upper_bound "upper_bound"(const T&)

cdef int cigar_exists(r)

cdef tuple clip_sizes(r)

cdef tuple clip_sizes_hard(r)

cdef int cigar_clip(r, int clip_length)

cpdef int is_overlapping(int x1, int x2, int y1, int y2) nogil


cdef float min_fractional_overlapping(int x1, int x2, int y1, int y2)

cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil

cdef bint span_position_distance(int x1, int x2, int y1, int y2, float norm, float thresh, ReadEnum_t read_enum, bint paired_end, int cigar_len1, int cigar_len2, bint trust_ins_len) nogil

cdef float position_distance(int x1, int x2, int y1, int y2) nogil

cdef class EventResult:
    """Data holder for classifying alignments into SV types"""
    cdef public int32_t contig_ref_start, contig_ref_end, contig2_ref_start, contig2_ref_end, contig_lc, contig_rc, contig2_lc, contig2_rc, \
        grp_id, event_id, n_expansion, stride, ref_poly_bases
    cdef public np.float32_t contig_left_weight, contig_right_weight, contig2_left_weight, contig2_right_weight, ref_rep, compress

    cdef public int32_t su, pe, supp, sc, NP, maxASsupp, plus, minus, spanning, double_clips, n_unmapped_mates, n_small_tlen, bnd, ras, fas, cipos95A, cipos95B
    cdef public np.float32_t NMpri, NMsupp, MAPQpri, MAPQsupp, NMbase, n_sa, n_xa, n_gaps

    cdef public int32_t posA, posB, svlen, query_gap, query_overlap, block_edge, ref_bases, remap_score, bad_clip_count, remap_ed, n_in_grp
    cdef public np.float32_t jitter, sqc, scw, clip_qual_ratio, outer_cn, inner_cn, fcc, rep, rep_sc, gc, neigh, neigh10kb, raw_reads_10kb, mcov, strand_binom_t
    cdef public bint preciseA, preciseB, linked, modified, remapped
    cdef public int8_t svlen_precise
    cdef public object contig, contig2, svtype, join_type, chrA, chrB, exp_seq, sample, type, \
        partners, GQ, GT, kind, ref_seq, variant_seq, left_ins_seq, right_ins_seq, site_info
