#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8

from libc.stdint cimport uint32_t


cdef float soft_clip_qual_corr(reads)


cdef struct WindowRate:
    float rate
    int index


cdef void window_rate(WindowRate *result, uint32_t cigar_l, uint32_t *cigar_p, int index, int window_size, bint reverse)
