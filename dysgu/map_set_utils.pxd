#cython: language_level=3

from libcpp.vector cimport vector as cpp_vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.utility cimport pair
from libcpp.string cimport string as cpp_string

import numpy as np
cimport numpy as np

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK


from libc.stdint cimport uint32_t, uint8_t, uint64_t, uint16_t, int32_t, int8_t


ctypedef cpp_vector[int] int_vec_t
ctypedef cpp_pair[int, int] get_val_result

# ctypedef Py_Int2IntVecMap[int, int_vec_t] node_dict_t
# ctypedef Py_IntVec2IntMap[int_vec_t, int] node_dict2_r_t

from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators

from libcpp cimport bool
# from pysam.libcalignedsegment cimport bam1_t

ctypedef enum ReadEnum_t:
    DISCORDANT = 0
    SPLIT = 1
    DELETION = 2
    INSERTION = 3
    BREAKEND = 4

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


cdef class Py_CoverageTrack:
    cdef public int max_coverage, current_chrom
    cdef public str outpath
    cdef public object cov_array #np.ndarray[np.int32, ndim=1] cov_array
    cdef public object infile
    # cdef void add(self, object a)
    # cdef void set_cov_array(self, int rname)
    # cdef void write_track(self)


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


cdef extern from "wrap_map_set2.h" nogil:
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


cdef int cigar_exists(r)


cdef tuple clip_sizes(r)


cdef tuple clip_sizes_hard(r)


cdef int cigar_clip(r, int clip_length)


cpdef int is_overlapping(int x1, int x2, int y1, int y2) nogil


cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil


cdef bint span_position_distance(int x1, int x2, int y1, int y2, float norm, float thresh, ReadEnum_t read_enum, bint paired_end) nogil


cdef float position_distance(int x1, int x2, int y1, int y2) nogil


cdef class EventResult:
    """Data holder for classifying alignments into SV types"""
    cdef public int32_t contig_ref_start, contig_ref_end, contig2_ref_start, contig2_ref_end, grp_id, event_id, n_expansion, stride, ref_poly_bases
    cdef public float contig_left_weight, contig_right_weight, contig2_left_weight, contig2_right_weight, ref_rep, compress

    cdef public int32_t su, pe, supp, sc, NP, maxASsupp, plus, minus, spanning, double_clips, n_unmapped_mates, n_small_tlen, bnd, ras, fas, cipos95A, cipos95B
    cdef public float DP, DApri, DN, NMpri, DAsupp, NMsupp, MAPQpri, MAPQsupp, NMbase, n_sa, n_xa, n_gaps

    cdef public int32_t posA, posB, svlen, query_gap, query_overlap, block_edge, ref_bases, remap_score, bad_clip_count, remap_ed, n_in_grp
    cdef public float jitter, sqc, scw, outer_cn, inner_cn, fcc, rep, rep_sc, gc, neigh, neigh10kb, raw_reads_10kb, mcov
    cdef public bint preciseA, preciseB, linked, modified, remapped
    cdef public int8_t svlen_precise
    cdef public object contig, contig2, svtype, join_type, chrA, chrB, exp_seq, sample, type, partners, GQ, SQ, GT, kind
