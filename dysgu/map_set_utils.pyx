#cython: language_level=3

from libcpp.vector cimport vector as cpp_vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.string cimport string as cpp_string

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK
from libc.stdint cimport uint32_t


# ctypedef cpp_vector[int] int_vec_t
# ctypedef cpp_pair[int, int] get_val_result

# ctypedef Py_Int2IntVecMap[int, int_vec_t] node_dict_t
# ctypedef Py_IntVec2IntMap[int_vec_t, int] node_dict2_r_t


# cdef extern from "wrap_map_set2.h":
#     cdef cppclass Int2IntMap:
#         Int2IntMap()
#         void insert(int, int)
#         void erase(int)
#         int has_key(int)
#         int get(int)
#         get_val_result get_value(int key)
#         int size()

cdef class Py_Int2IntMap:
    """Fast 32bit integer to 32bit integer unordered map using tsl::robin-map"""
    # cdef Int2IntMap *thisptr
    def __cinit__(self):
        self.thisptr = new Int2IntMap()
    def __dealloc__(self):
        del self.thisptr
    cpdef void insert(self, int key, int value):
        self.thisptr.insert(key, value)
    cpdef void erase(self, int key):
        self.thisptr.erase(key)
    cpdef int has_key(self, int key):
        return self.thisptr.has_key(key)
    cpdef int get(self, int key):
        return self.thisptr.get(key)
    cpdef get_val_result get_value(self, int key):
        return self.thisptr.get_value(key)
    cpdef int size(self):
        return self.thisptr.size()


# cdef extern from "wrap_map_set2.h":
#     cdef cppclass IntSet:
#         IntSet()
#         void insert(int)
#         void erase(int)
#         int has_key(int)
#         int get(int)
#         int size()

cdef class Py_IntSet:
    """Fast 32 bit int set using tsl::robin-set"""
    # cdef IntSet *thisptr
    def __cinit__(self):
        self.thisptr = new IntSet()
    def __dealloc__(self):
        del self.thisptr
    cpdef void insert(self, int key):
        self.thisptr.insert(key)
    cpdef void erase(self, int key):
        self.thisptr.erase(key)
    cpdef int has_key(self, int key):
        return self.thisptr.has_key(key)
    cpdef int size(self):
        return self.thisptr.size()


cdef class Py_StrSet:

    def __cinit__(self):
        self.thisptr = new StrSet()
    def __dealloc__(self):
        del self.thisptr
    cpdef void insert(self, cpp_string key):
        self.thisptr.insert(key)
    cpdef void erase(self, cpp_string key):
        self.thisptr.erase(key)
    cpdef int has_key(self, cpp_string key):
        return self.thisptr.has_key(key)
    cpdef int size(self):
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


cdef int cigar_clip(r, int clip_length):

    c = r.cigartuples
    if not c:
        return 0
    if (c[0][0] == 4 and c[0][1] > clip_length) or (c[-1][0] == 4 and c[-1][1] > clip_length):
        return 1
    return 0

