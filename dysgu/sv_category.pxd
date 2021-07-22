#cython: language_level=3

cdef class AlignmentItem:
    """Data holder for classifying alignments into SV types"""
    cdef public int chrA, chrB, priA, priB, rA, rB, posA, endA, posB, endB, strandA, strandB, left_clipA, right_clipA,\
        left_clipB, right_clipB, breakA_precise, breakB_precise, breakA, breakB, a_qstart, a_qend, b_qstart, b_qend,\
        a_len, b_len, query_gap, read_overlaps_mate, size_inferred, query_overlap, inferred_sv_len
    cdef public str svtype, join_type
    cdef public object read_a, read_b


cdef void classify_d(AlignmentItem v)