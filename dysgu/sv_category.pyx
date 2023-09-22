#cython: language_level=3

from dysgu.map_set_utils cimport is_overlapping
from dysgu.map_set_utils import echo


cdef class AlignmentItem:
    """Data holder for classifying alignments into SV types"""
    # cdef public int chrA, chrB, priA, priB, rA, rB, posA, endA, posB, endB, strandA, strandB, left_clipA, right_clipA,\
    #     left_clipB, right_clipB, breakA_precise, breakB_precise, breakA, breakB, a_qstart, a_qend, b_qstart, b_qend,\
    #     a_len, b_len, query_gap, read_overlaps_mate, size_inferred, query_overlap, inferred_sv_len
    # cdef public str svtype, join_type
    # cdef public object read_a, read_b
    def __cinit__(self, int chrA, int chrB, int priA, int priB, int rA, int rB, int posA, int endA, int posB, int endB,
                  int strandA, int strandB, int left_clipA, int right_clipA, int left_clipB, int right_clipB,
                  int a_qstart, int a_qend, int b_qstart, int b_qend, int a_len, int b_len,
                  int read_overlaps_mate, object read_a,
                  object read_b):
        self.chrA = chrA
        self.chrB = chrB
        self.priA = priA
        self.priB = priB
        self.rA = rA
        self.rB = rB
        self.posA = posA
        self.endA = endA
        self.posB = posB
        self.endB = endB
        self.strandA = strandA
        self.strandB = strandB
        self.left_clipA = left_clipA
        self.right_clipA = right_clipA
        self.left_clipB = left_clipB
        self.right_clipB = right_clipB
        self.a_qstart = a_qstart
        self.a_qend = a_qend
        self.b_qstart = b_qstart
        self.b_qend = b_qend
        self.a_len = a_len
        self.b_len = b_len
        self.read_overlaps_mate = read_overlaps_mate
        self.read_a = read_a
        self.read_b = read_b

        # These will be determined
        self.breakA_precise = 0
        self.breakB_precise = 0
        self.breakA = -1
        self.breakB = -1
        self.svtype = ""
        self.join_type = ""
        self.query_gap = 0
        self.inferred_sv_len = -1
        self.size_inferred = 0  # set to 1 if insertion size was inferred


cdef void two_primary(AlignmentItem v):

    if v.posA < v.posB or (v.posA == v.posB and v.endA < v.endB):

        if v.strandA == 3 and v.strandB == 5:  # DEL type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            if v.endA >= v.posB:
                v.svtype = "INS"
                v.join_type = "3to5"

            else:
                v.svtype = "DEL"
                v.join_type = "3to5"

        elif v.strandA == 5 and v.strandB == 3:  # DUP type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

        else:  # INV type
            # Break to left or right
            if v.strandA == 5:

                if not (v.left_clipA or v.left_clipB) and (v.right_clipA or v.right_clipB):
                    v.breakA = v.endA
                    if v.right_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.endB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "3to3"

                else:
                    v.breakA = v.posA
                    if v.left_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.posB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "5to5"

            else:  # strand == 3
                if not (v.right_clipA or v.right_clipB) and (v.left_clipA or v.left_clipB):  # only left clips

                    if (v.posA < v.posB and v.left_clipB) or (v.posB < v.posA and v.left_clipA):  # guess right
                        v.breakA = v.endA
                        if v.right_clipA:
                            v.breakA_precise = 1

                        v.breakB = v.endB
                        if v.right_clipB:
                            v.breakB_precise = 1

                        v.svtype = "INV"
                        v.join_type = "3to3"

                    else:  # guess left
                        v.breakA = v.posA
                        if v.left_clipA:
                            v.breakA_precise = 1

                        v.breakB = v.posB
                        if v.left_clipB:
                            v.breakB_precise = 1

                        v.svtype = "INV"
                        v.join_type = "5to5"

                else:  # either no clips or right clips

                    v.breakA = v.endA
                    if v.left_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.endB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "3to3"

    else:  # B < A

        if v.strandA == 5 and v.strandB == 3:  # DEL type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            if v.endB >= v.posA:
                v.svtype = "INS"
                v.join_type = "3to5"

            else:
                v.svtype = "DEL"
                v.join_type = "3to5"

        elif v.strandB == 5 and v.strandA == 3:  # DUP type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

            # INV type
        elif v.strandA == 5:

            if not (v.left_clipA or v.left_clipB) and (v.right_clipA or v.right_clipB):
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "3to3"

            else:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.posB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "5to5"

        elif v.strandA == 3:

            if not (v.right_clipA or v.right_clipB) and (v.left_clipA or v.left_clipB):
                v.breakA = v.posA
                if v.right_clipA:
                    v.breakA_precise = 1

                v.breakB = v.posB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "5to5"

            else:

                v.breakA = v.endA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "3to3"


cdef void same_read(AlignmentItem v):

    cdef int query_gap = 0
    cdef int ref_gap = 0
    cdef int query_overlap = 0
    cdef int inferred_sv_len = -1

    if v.posA < v.posB or (v.posA == v.posB and v.endA <= v.endB):  # A is first

        if v.strandA == v.strandB:  # Same for 3 and 5 strand reads

            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested DUP
                if v.left_clipA:
                    v.breakA = v.posA
                    v.breakA_precise = 1
                elif v.right_clipA:
                    v.breakA = v.endA
                    v.breakA_precise = 1
                else:
                    v.breakA = v.endA

                if v.left_clipB:
                    v.breakB = v.posB
                    v.breakB_precise = 1
                elif v.right_clipB:
                    v.breakB = v.endB
                    v.breakB_precise = 1
                else:
                    v.breakB = v.endB

                v.svtype = "DUP"
                # v.svtype = "INS"
                # Check if gap on query is bigger than gap on reference
                ref_gap = abs(v.breakB - v.breakA)
                if v.a_qstart > v.b_qstart:  # B is first on query seq
                    align_overlap = v.b_qend - v.a_qstart  # note, this is negative if there is a gap between alignments
                else:
                    align_overlap = v.a_qend - v.b_qstart
                # align_overlap = 0 if align_overlap < 0 else align_overlap
                v.inferred_sv_len = ref_gap - align_overlap
                if align_overlap < 0:
                    v.query_gap = align_overlap
                else:
                    v.query_overlap = align_overlap
                v.join_type = "5to3"

            elif not v.left_clipA and not v.right_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1

                # Check if gap on query is bigger than gap on reference; call insertion if it is
                query_gap = v.b_qstart - v.a_qend  # as A is first
                ref_gap = abs(v.breakB - v.breakA)
                if query_gap < 0:
                    ref_gap += abs(query_gap)
                    v.query_overlap = abs(query_gap)
                    v.query_gap = 0

                else:
                    v.query_overlap = 0
                    v.query_gap = query_gap

                if ref_gap < query_gap:
                    v.svtype = "INS"
                    v.join_type = "3to5"
                    v.inferred_sv_len = query_gap

                elif v.read_overlaps_mate:
                    v.svtype = "DUP"
                    v.join_type = "5to3"
                else:
                    v.svtype = "DEL"
                    v.join_type = "3to5"
                    v.inferred_sv_len = ref_gap

            elif not v.right_clipA and not v.left_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif not v.left_clipA and not v.left_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            else:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"

        elif v.strandA == 5 and v.strandB == 3:

            # Break right
            if not v.left_clipA and not v.left_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "3to3"

            # Break left
            elif not v.right_clipA and not v.right_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "5to5"

            elif is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Inverted duplication
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif v.right_clipA and v.left_clipB:
                v.breakA = v.endA
                v.breakB = v.posB
                v.svtype = "INV"
                v.join_type = "3to5"

            else:
            # elif v.left_clipA and v.right_clipB:
                v.breakA = v.posA
                v.breakB = v.endB
                v.svtype = "INV"
                v.join_type = "5to3"

        else:  # INV type
            # Break left
            if v.left_clipA and v.left_clipB:
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                v.join_type = "5to5"
            # Break right
            elif v.right_clipA and v.right_clipB:
                v.breakA = v.endA
                v.breakB = v.endB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            else:  # Guess using pos only
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                if v.strandA == 5:
                    v.join_type = "5to5"
                else:
                    v.join_type = "3to3"

    # B is first or B == A
    else:

        if v.strandA == v.strandB:

            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested DUP
                if v.left_clipA:
                    v.breakA = v.posA
                    v.breakA_precise = 1
                elif v.right_clipA:
                    v.breakA = v.endA
                    v.breakA_precise = 1
                else:
                    v.breakA = v.endA

                if v.left_clipB:
                    v.breakB = v.posB
                    v.breakB_precise = 1
                elif v.right_clipB:
                    v.breakB = v.endB
                    v.breakB_precise = 1
                else:
                    v.breakB = v.endB

                v.svtype = "INS"
                # Check if gap on query is bigger than gap on reference
                ref_gap = abs(v.breakA - v.breakB)

                if v.a_qstart > v.b_qstart:  # B is first on query seq
                    align_overlap = v.b_qend - v.a_qstart  # note, this is negative if there is a gap between alignments
                else:
                    align_overlap = v.a_qend - v.b_qstart

                v.inferred_sv_len = ref_gap - align_overlap

                if align_overlap < 0:
                    v.query_gap = align_overlap
                else:
                    v.query_overlap = align_overlap

                v.join_type = "5to3"

            elif not v.left_clipB and not v.right_clipA:
                v.breakA = v.posA
                v.breakB = v.endB
                v.breakA_precise = 1
                v.breakB_precise = 1

                # Check if gap on query is bigger than gap on reference; call insertion if it is
                query_gap = v.b_qend - v.a_qstart  # as B is first
                ref_gap = abs(v.breakA - v.breakB)

                if query_gap < 0:
                    ref_gap += abs(query_gap)
                    v.query_overlap = abs(query_gap)
                    v.query_gap = 0

                else:
                    v.query_overlap = 0
                    v.query_gap = query_gap

                if ref_gap < query_gap:
                    v.svtype = "INS"
                    v.join_type = "3to5"
                    v.inferred_sv_len = query_gap

                elif v.read_overlaps_mate:
                    v.svtype = "DUP"
                    v.join_type = "5to3"
                else:
                    v.svtype = "DEL"
                    v.join_type = "3to5"
                    v.inferred_sv_len = ref_gap

            elif not v.right_clipB and not v.left_clipA:
                v.breakA = v.endA
                v.breakB = v.posB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            else:
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "BND"
                v.join_type = f"{v.strandA}to{v.strandB}"

        else:  # INV type
            # Break left
            if v.left_clipA and v.left_clipB:
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                v.join_type = "5to5"

            # Break right
            elif v.right_clipA and v.right_clipB:
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.breakA = v.endA
                v.breakB = v.endB
                v.svtype = "INV"
                v.join_type = "3to3"

            else:  # Guess using pos only
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                if v.strandA == 5:
                    v.join_type = "5to5"
                else:
                    v.join_type = "3to3"


cdef void different_read(AlignmentItem v):

    if v.posA < v.posB or (v.posA == v.posB and v.endA <= v.endB):  # A is first

        if v.strandA == 3 and v.strandB == 5:  # DEL type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "INS"
            v.join_type = "3to5"

        elif v.strandA == 5 and v.strandB == 3:  # DUP type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

        elif v.strandA == v.strandB:
            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested

                if v.strandA == 3:  # Both forward strand
                    v.breakA = v.endA
                    if v.right_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.endB
                    if v.right_clipB:
                        v.breakB_precise = 1
                    v.svtype = "INV"
                    v.join_type = "3to3"

                else:  # Both forward strand
                    v.breakA = v.posA
                    if v.left_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.posB
                    if v.left_clipB:
                        v.breakB_precise = 1
                    v.svtype = "INV"
                    v.join_type = "5to5"

            elif not v.right_clipA and not v.right_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"

            elif not v.left_clipA and not v.left_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            elif v.right_clipA and v.left_clipB:
                v.breakA = v.endA
                v.breakA_precise = 1
                v.breakB = v.posB
                v.breakB_precise = 1
                # v.svtype = "INV:DUP"
                v.svtype = "INV"
                v.join_type = "5to3"

            else:
                v.breakA = v.posA
                v.breakA_precise = 1
                v.breakB = v.endB
                v.breakB_precise = 1
                v.svtype = "INV"
                # v.svtype = "INV:DUP"
                v.join_type = "3to5"

    else:  # B is first; B <= A

        if v.strandA == 5 and v.strandB == 3:  # DEL type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "INS"
            v.join_type = "3to5"

        elif v.strandA == 3 and v.strandB == 5:  # DUP type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

        elif v.strandA == v.strandB:  # INV type

            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested DUP
                if v.left_clipA:
                    v.breakA = v.posA
                    v.breakA_precise = 1
                elif v.right_clipA:
                    v.breakA = v.endA
                    v.breakA_precise = 1
                else:
                    v.breakA = v.endA

                if v.left_clipB:
                    v.breakB = v.posB
                    v.breakB_precise = 1
                elif v.right_clipB:
                    v.breakB = v.endB
                    v.breakB_precise = 1
                else:
                    v.breakB = v.endB
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif v.right_clipA and v.left_clipB:
                v.breakA = v.endA
                v.breakB = v.posB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DEL"
                v.join_type = "3to5"

            elif v.left_clipA and v.right_clipB:
                v.breakA = v.posA
                v.breakB = v.endB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif not v.left_clipB and not v.left_clipA:

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            elif not v.right_clipB and not v.right_clipA:

                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1

                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"

            elif v.strandA == 3:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            else:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"


cdef void translocation(AlignmentItem v):

    v.svtype = "TRA"

    if v.left_clipA:
        v.breakA = v.posA
        v.breakA_precise = 1
    elif v.right_clipA:
        v.breakA = v.endA
        v.breakA_precise = 1
    else:
        v.breakA = v.posA
        v.breakA_precise = 1

    if v.left_clipB:
        v.breakB = v.posB
        v.breakB_precise = 1
    elif v.right_clipB:
        v.breakB = v.endB
        v.breakB_precise = 1
    else:
        v.breakB = v.posB
        v.breakB_precise = 1
    v.join_type = f"{v.strandA}to{v.strandB}"

    cdef int query_gap = 0
    cdef int query_overlap = 0
    if v.rA == v.rB:  # same read
        if v.b_qstart < v.a_qstart:  # B is first on query
            query_gap = v.a_qstart - v.b_qend
        else:
            query_gap = v.b_qstart - v.a_qend
        if query_gap < 0:
            query_overlap = abs(query_gap)
            query_gap = 0

    v.query_gap = query_gap
    v.query_overlap = query_overlap


cdef void classify_d(AlignmentItem v):  #, debug=False):

    v.breakA_precise = 0
    v.breakB_precise = 0
    if v.chrA != v.chrB:
        translocation(v)
    else:  # Intra-chromosomal. Find join type first, different for pri-sup, pri-pri
        if v.priA and v.priB:  # Both primary
            two_primary(v)
        else:  # One is a supplementary
            if v.rA == v.rB:
                same_read(v)
            else:
                different_read(v)


    # Debug
    # if debug:
    #     echo(v.read_a.qname, v.rA == v.rB, v.priA and v.priB, v.inferred_sv_len, "strands", v.strandA == v.strandB, "pos", v.posA, v.posB, "breaks", v.breakA, v.breakB, "a_q > b_q", v.a_qstart > v.b_qstart)
