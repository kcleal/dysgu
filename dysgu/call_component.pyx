#cython: language_level=3, boundscheck=True, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

from __future__ import absolute_import
from collections import Counter, defaultdict
import click
import numpy as np
cimport numpy as np
from numpy.random import normal
import itertools
from dysgu import assembler
from dysgu.map_set_utils import echo
from dysgu.map_set_utils cimport hash as xxhasher
from dysgu.map_set_utils cimport is_overlapping, clip_sizes_hard, EventResult
from dysgu.post_call_metrics cimport soft_clip_qual_corr

import warnings
from pysam.libcalignedsegment cimport AlignedSegment
from libc.stdint cimport uint64_t

from pysam.libchtslib cimport bam_get_qname

from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster


warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

np.random.seed(1)


ctypedef EventResult EventResult_t


cdef class AlignmentItem:
    """Data holder for classifying alignments into SV types"""
    cdef public int chrA, chrB, priA, priB, rA, rB, posA, endA, posB, endB, strandA, strandB, left_clipA, right_clipA,\
        left_clipB, right_clipB, breakA_precise, breakB_precise, breakA, breakB, a_qstart, a_qend, b_qstart, b_qend,\
        query_gap, read_overlaps_mate, size_inferred, query_overlap, inferred_sv_len
    cdef public str svtype, join_type
    cdef public object read_a, read_b
    def __cinit__(self, int chrA, int chrB, int priA, int priB, int rA, int rB, int posA, int endA, int posB, int endB,
                  int strandA, int strandB, int left_clipA, int right_clipA, int left_clipB, int right_clipB,
                  int a_qstart, int a_qend, int b_qstart, int b_qend, int read_overlaps_mate, object read_a,
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


cdef n_aligned_bases(ct):
    cdef int opp, l, aligned, large_gaps, n_small_gaps
    aligned, large_gaps, n_small_gaps = 0, 0, 0
    for opp, l in ct:
        if opp == 0:
            aligned += l
        elif opp == 1 or opp == 2:
            if l >= 30:
                large_gaps += l
            else:
                n_small_gaps += 1

    return float(aligned), float(large_gaps), float(n_small_gaps)


cdef count_attributes2(reads1, reads2, spanning, int extended_tags, float insert_ppf, generic_ins,
                       EventResult_t er):

    cdef float DP = 0
    cdef float DApri = 0
    cdef float DN = 0
    cdef float DAsupp = 0

    cdef float NMpri = 0
    cdef float NMsupp = 0
    cdef int maxASsupp = 0
    cdef float MAPQpri = 0
    cdef float MAPQsupp = 0
    cdef float NMbase = 0
    cdef int n_sa = 0
    cdef int n_xa = 0
    cdef float n_gaps = 0

    cdef int total_pri = 0

    paired_end = set([])
    seen = set([])

    er.spanning = len(spanning)
    er.bnd = len(generic_ins)

    cdef int flag, index, this_as
    cdef float a_bases, large_gaps, n_small_gaps
    # extended tags: "DP", "DN", "DApri", "DAsupp",
    for index, a in enumerate(itertools.chain(reads1, reads2)):

        qname = a.qname
        if qname not in seen:
            seen.add(qname)

        if extended_tags:
            if a.has_tag("ZP"):
                DP += float(a.get_tag("ZP"))
            if a.has_tag("ZN"):
                DN += float(a.get_tag("ZN"))

        flag = a.flag
        if flag & 2:
            er.NP += 1
        if a.flag & 1 and a.tlen and abs(a.tlen) < insert_ppf:
            er.n_small_tlen += 1
        if paired_end and flag & 8:
            er.n_unmapped_mates += 1

        left_clip, right_clip = clip_sizes_hard(a)
        if left_clip > 0 and right_clip > 0:
            er.double_clips += 1

        has_sa = a.has_tag("SA")
        if has_sa:
            n_sa += a.get_tag("SA").count(";")
        if a.has_tag("XA"):
            n_xa += a.get_tag("XA").count(";")

        a_bases, large_gaps, n_small_gaps = n_aligned_bases(a.cigartuples)
        if a_bases > 0:
            n_gaps += n_small_gaps / a_bases

        if flag & 2304:  # Supplementary (and not primary if -M if flagged using bwa)
            er.supp += 1
            MAPQsupp += a.mapq
            if extended_tags and a.has_tag("ZA"):
                DAsupp += float(a.get_tag("ZA"))
            if a.has_tag("NM"):
                if a_bases:
                    NMsupp += float(a.get_tag("NM")) / a_bases

            if a.has_tag("AS"):
                this_as = a.get_tag("AS")
                if this_as > maxASsupp:
                    maxASsupp = this_as

        else:  # Primary reads
            total_pri += 1

            MAPQpri += a.mapq
            if index >= len(reads2) and qname in paired_end:  # If two primary reads from same pair
                er.pe += 1
            else:
                paired_end.add(qname)

            if extended_tags and a.has_tag("ZA"):
                DApri += float(a.get_tag("ZA"))
            if a.has_tag("NM"):
                if a_bases:
                    nm = float(a.get_tag("NM"))
                    NMpri += nm / a_bases
                    NMbase += (nm - large_gaps) / a_bases

        if flag & 16:
            er.minus += 1
        else:
            er.plus += 1

        ct = a.cigartuples
        if ct[0][0] == 4 or ct[-1][0] == 4:
            er.sc += 1

    for a in spanning:

        qname = a.qname
        if qname not in seen:
            seen.add(qname)

        if extended_tags:
            if a.has_tag("ZP"):
                DP += float(a.get_tag("ZP"))
            if a.has_tag("ZN"):
                DN += float(a.get_tag("ZN"))

        flag = a.flag
        if flag & 2:
            er.NP += 1

        left_clip, right_clip = clip_sizes_hard(a)
        if left_clip > 0 and right_clip > 0:
            er.double_clips += 1

        if a.has_tag("SA"):
            n_sa += a.get_tag("SA").count(";")
        if a.has_tag("XA"):
            n_xa += a.get_tag("XA").count(";")

        a_bases, large_gaps, n_small_gaps = n_aligned_bases(a.cigartuples)
        if a_bases > 0:
            n_gaps += n_small_gaps / a_bases

        if flag & 2304:  # Supplementary
            er.supp += 1
            MAPQsupp += a.mapq
            if extended_tags and a.has_tag("ZA"):
                DAsupp += float(a.get_tag("ZA"))
            if a.has_tag("NM"):
                if a_bases:
                    NMsupp += float(a.get_tag("NM")) / a_bases

            if a.has_tag("AS"):
                this_as = a.get_tag("AS")
                if this_as > maxASsupp:
                    maxASsupp = this_as

        else:  # Primary reads
            MAPQpri += a.mapq
            total_pri += 1
            if qname in paired_end:  # If two primary reads from same pair
                er.pe += 2
            else:
                paired_end.add(qname)
            if extended_tags and a.has_tag("ZA"):
                DApri += float(a.get_tag("ZA"))
            if a.has_tag("NM"):
                if a_bases:
                    nm = float(a.get_tag("NM"))
                    NMpri += nm / a_bases
                    NMbase += (nm - large_gaps) / a_bases

        if flag & 16:
            er.minus += 1
        else:
            er.plus += 1

    cdef int tot = er.supp + total_pri

    er.NMpri = (NMpri / total_pri) * 100 if total_pri > 0 else 0
    er.NMsupp = (NMsupp / er.supp) * 100 if er.supp > 0 else 0
    er.NMbase = (NMbase / tot) * 100 if tot > 0 else 0
    er.n_gaps = (n_gaps / tot) * 100 if tot > 0 else 0

    if extended_tags:
        er.DP = DP / tot if tot > 0 else 0
        er.DN = DN / tot if tot > 0 else 0
        er.DApri = DApri / total_pri if total_pri > 0 else 0
        er.DAsupp = DAsupp / er.supp if er.supp > 0 else 0

    er.MAPQpri = MAPQpri / total_pri if total_pri > 0 else 0
    er.MAPQsupp = MAPQsupp / er.supp if er.supp > 0 else 0
    er.n_sa = n_sa / tot if tot > 0 else 0
    er.n_xa = n_xa / tot if tot > 0 else 0

    er.maxASsupp = maxASsupp

    if len(reads2) == 0:
        er.pe = len(reads1)
    er.su = er.pe + er.supp + (2*er.spanning) + er.bnd


cdef int within_read_end_position(event_pos, svtype, cigartuples, cigar_index):
    # If a deletion has a poorly aligned section at the end of pos2, merge nearby events until target matches
    cdef int i, end, idx, opp, length, cigar_skip, n_aligned_bases, target_bases
    if svtype == 1:  # insertion (or other e.g. duplication/inversion within read)
        return event_pos + 1
    else:  # deletion type
        end = event_pos + cigartuples[cigar_index][1]
        return end


cdef guess_informative_pair(aligns):

    if len(aligns) == 2:
        a_cigar_info, a = aligns[0]
        b_cigar_info, b = aligns[1]

        # check for paired-end read through with no SA tag
        if a.flag & 1 and a_cigar_info[5] == -1 and b_cigar_info[5] == -1:
            if a.pos == b.pos and a.reference_end == b.reference_end:
                extent_left_same = True
                extent_right_same = True
                if a.cigartuples[0][0] == 4 and not (a.cigartuples[0] == b.cigartuples[0]):
                    extent_left_same = False
                if extent_left_same and a.cigartuples[-1][0] == 4 and not (a.cigartuples[-1] == b.cigartuples[-1]):
                    extent_right_same = False
                if extent_left_same and extent_right_same:
                    return None

        # within read sv
        if 0 < a_cigar_info[5] < len(a.cigartuples) - 1:
            cigar_index = a_cigar_info[5]
            event_pos = a_cigar_info[6]
            ci = a.cigartuples[cigar_index]
            return (ci[0],
                    a.rname,
                    event_pos,
                    event_pos + 1 if ci[0] == 1 else event_pos + a.cigartuples[cigar_index][1] + 1,
                    a.cigartuples[cigar_index][1],
                    a)

        elif 0 < b_cigar_info[5] < len(b.cigartuples) - 1:
            cigar_index = b_cigar_info[5]
            event_pos = b_cigar_info[6]
            ci = b.cigartuples[cigar_index]
            return (ci[0],
                    b.rname,
                    event_pos,
                    event_pos + 1 if ci[0] == 1 else event_pos + b.cigartuples[cigar_index][1] + 1,
                    b.cigartuples[cigar_index][1],
                    b)

        a_pos = a_cigar_info[6]  # Position may have been inferred from SA tag, use this if available
        b_pos = b_cigar_info[6]
        if a_pos == b_pos:
            # make sure different breaks are mapped
            if (a.cigartuples[0][0] == 4 and b.cigartuples[0][0] == 4) or (a.cigartuples[-1][0] == 4 and b.cigartuples[-1][0] == 4):
                return None

        # a and b will be same on same chrom
        if a_pos < b_pos:
            return a, b

        return b, a

    pri_first = None
    sup_first = None
    pri_second = None
    sup_second = None
    for node_info, i in aligns:

        if not i.flag & 2304:  # (Not pri, supplementary) --> is primary
            if i.flag & 64:
                pri_first = i
            else:
                pri_second = i
        else: # Supplementary, -M flag of bwa marks supplementary as not primary
            if i.flag & 64:
                sup_first = i
            else:
                sup_second = i
    a = None
    b = None
    if pri_first and sup_first:
        a = pri_first
        b = sup_first
    elif pri_second and sup_second:
        a = pri_second
        b = sup_second
    elif pri_first and pri_second:
        a = pri_first
        b = pri_second
    if a is None:
        return None
    if a.pos < b.pos:
        return a, b
    return b, a


cdef int same_read_overlaps_mate(a_chrom, b_chrom, a_start, a_end, b_start, b_end, a, b):
    # If data is paired-end, check if one supplementary overlaps the other primary read
    aflag = a.flag
    bflag = b.flag
    if aflag & 1 and not aflag & 8 and a_chrom == b_chrom and aflag & 64 == bflag & 64:  # same read, one is supplementary
        if is_overlapping(b_start, b_end, a_start, a_end):
            return 1
    return 0


cdef make_sv_alignment_item(a, b):

    a_ct = a.cigartuples
    b_ct = b.cigartuples

    a_qstart, a_qend, b_qstart, b_qend = start_end_query_pair(a, b)

    # Soft-clips for the chosen pair, plus template start of alignment
    left_clip_a, right_clip_a, left_clip_b, right_clip_b = mask_soft_clips(a.flag, b.flag, a_ct, b_ct)

    a_chrom, b_chrom = a.rname, b.rname
    a_start, a_end = a.pos, a.reference_end
    b_start, b_end = b.pos, b.reference_end
    read_overlaps_mate = same_read_overlaps_mate(a_chrom, b_chrom, a_start, a_end, b_start, b_end, a, b)

    v_item = AlignmentItem(a_chrom, b_chrom,
                           int(not a.flag & 2304),  # is primary
                           int(not b.flag & 2304),
                           1 if a.flag & 64 else 2,
                           1 if b.flag & 64 else 2,
                           a_start, a_end,
                           b_start, b_end,
                           3 if a.flag & 16 == 0 else 5,
                           3 if b.flag & 16 == 0 else 5,
                           left_clip_a, right_clip_a,
                           left_clip_b, right_clip_b,
                           a_qstart, a_qend, b_qstart, b_qend,
                           read_overlaps_mate,
                           a, b
                           )

    if v_item.left_clipA and v_item.right_clipA:
        if a_ct[0][1] >= a_ct[-1][1]:
            v_item.right_clipA = 0
        else:
            v_item.left_clipA = 0

    if v_item.left_clipB and v_item.right_clipB:
        if b_ct[0][1] >= b_ct[-1][1]:
            v_item.right_clipB = 0
        else:
            v_item.left_clipB = 0

    return v_item


cdef make_generic_insertion_item(aln, int insert_size, int insert_std):

    if aln.flag & 2304:  # skip supplementary
        return None
    v_item = make_sv_alignment_item(aln, aln)
    v_item.read_b = None
    v_item.read_overlaps_mate = 0

    cdef int dist_to_break = 0
    cdef int rand_insert_pos = 0
    if v_item.left_clipA:  # use clip to guess break point
        v_item.breakA = aln.pos
        v_item.breakB = aln.pos
        v_item.breakA_precise = 1
        v_item.breakB_precise = 1
        v_item.join_type = "5to3"

    elif v_item.right_clipA:
        v_item.breakA = aln.reference_end
        v_item.breakB = aln.reference_end
        v_item.breakA_precise = 1
        v_item.breakB_precise = 1
        v_item.join_type = "3to5"

    else:
        if aln.flag & 1 and aln.flag & 2: # paired and normal primary pairing
            if aln.flag & 16:  # guess left
                v_item.breakA = aln.pos
                v_item.breakB = aln.pos
                v_item.join_type = "5to3"
            else:
                v_item.breakA = aln.reference_end
                v_item.breakB = aln.reference_end
                v_item.join_type = "3to5"
        else:
            return None

    v_item.svtype = "INS"

    # Try and guess the insertion size using insert size and distance to break
    dist_to_break = aln.reference_end - aln.pos
    v_item.size_inferred = 1

    if insert_std > 0:
        rand_insert_pos = insert_size - dist_to_break + int(normal(0, insert_std))
    else:  # single read mode
        rand_insert_pos = 100
    v_item.inferred_sv_len = 0 if rand_insert_pos < 0 else rand_insert_pos

    return v_item


cpdef int assign_contig_to_break(asmb, EventResult_t er, side, spanning):
    if not asmb:
        return 0
    cdef int ref_bases = 0
    if spanning:
        er.contig = asmb["contig"]
        ref_bases += asmb["ref_bases"]
        er.contig_ref_start = asmb["ref_start"]
        er.contig_ref_end = asmb["ref_end"]
        er.contig_left_weight = 0
        er.contig_right_weight = 0
    elif asmb["left_clips"] or asmb["right_clips"]:
        if asmb["left_clips"] > asmb["right_clips"]:
            asmb_pos = asmb["ref_start"]
        else:
            asmb_pos = asmb["ref_end"]
        if side == "A":
            other_cont = "contig2"
            other_pos = er.posB
            other_chrom = er.chrB
            current_cont = "contig"
            current_pos = er.posA
            current_chrom = er.chrA
        else:
            other_cont = "contig"
            other_pos = er.posA
            other_chrom = er.chrA
            current_cont = "contig2"
            current_pos = er.posB
            current_chrom = er.chrB
        if other_chrom == current_chrom and abs(asmb_pos - other_pos) < abs(asmb_pos - current_pos):
            # Assign contig to opposite side
            if side == "A":
                er.contig2 = asmb["contig"]
                er.contig2_ref_start = asmb["ref_start"]
                er.contig2_ref_end = asmb["ref_end"]
                er.contig2_left_weight = asmb["left_weight"]
                er.contig2_right_weight = asmb["right_weight"]
            else:
                er.contig = asmb["contig"]
                er.contig_ref_start = asmb["ref_start"]
                er.contig_ref_end = asmb["ref_end"]
                er.contig_left_weight = asmb["left_weight"]
                er.contig_right_weight = asmb["right_weight"]

        else:
            # assign contigs to current side
            if side == "A":
                er.contig = asmb["contig"]
                er.contig_ref_start = asmb["ref_start"]
                er.contig_ref_end = asmb["ref_end"]
                er.contig_left_weight = asmb["left_weight"]
                er.contig_right_weight = asmb["right_weight"]
            else:
                er.contig2 = asmb["contig"]
                er.contig2_ref_start = asmb["ref_start"]
                er.contig2_ref_end = asmb["ref_end"]
                er.contig2_left_weight = asmb["left_weight"]
                er.contig2_right_weight = asmb["right_weight"]

        ref_bases += asmb["ref_bases"]
    return ref_bases


cdef make_single_call(sub_informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble, spanning_alignments,
                      extended_tags, svlen_precise, generic_ins):

    cdef EventResult_t er = EventResult()

    precise_a = []
    precise_b = []
    u_reads = []
    v_reads = []

    if len(sub_informative) == 0:
        si = generic_ins
    else:
        si = sub_informative

    for v_item in si:
        if v_item.breakA_precise:
            precise_a.append(v_item.breakA)
        if v_item.breakB_precise:
            precise_b.append(v_item.breakB)
        u_reads.append(v_item.read_a)
        if v_item.read_b is not None:
            v_reads.append(v_item.read_b)

    call_informative = Counter([(itm.svtype, itm.join_type) for itm in si]).most_common()
    svtype, jointype = call_informative[0][0]

    make_call(si, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev, er)

    min_found_support = call_informative[0][1]

    if len(sub_informative) > 0:
        count_attributes2(u_reads, v_reads, [], extended_tags, insert_ppf, generic_ins, er)
    else:
        count_attributes2([], [], [], extended_tags, insert_ppf, generic_ins, er)

    er.contig = None
    er.contig_left_weight = 0
    er.contig_right_weight = 0
    er.contig2 = None
    er.contig2_left_weight = 0
    er.contig2_right_weight = 0
    er.contig_ref_start = -1
    er.contig_ref_end = -1
    er.contig2_ref_start = -1
    er.contig2_ref_end = -1
    as1 = None
    as2 = None
    ref_bases = 0
    if to_assemble or len(spanning_alignments) > 0:
        if er.preciseA:
            as1 = assembler.base_assemble(u_reads, er.posA, 500)
            ref_bases += assign_contig_to_break(as1, er, "A", spanning_alignments)
        if er.preciseB:
            as2 = assembler.base_assemble(v_reads, er.posB, 500)
            ref_bases += assign_contig_to_break(as2, er, "B", 0)
    er.linked = 0
    er.block_edge = 0
    er.ref_bases = ref_bases
    er.svlen_precise = svlen_precise  # if 0 then soft-clip will be remapped
    if er.svlen_precise == 0:
        corr_score = soft_clip_qual_corr(u_reads + v_reads)
        er.sqc = corr_score
    else:
        er.sqc = -1

    return er


cdef partition_single(informative, insert_size, insert_stdev, insert_ppf, spanning_alignments,
                      min_support, extended_tags, to_assemble, generic_insertions):

    # spanning alignments is empty
    cdef AlignmentItem v_item
    cdef int idx = 0
    cdef np.ndarray[double, ndim=2] coords = np.zeros((len(informative), 2))
    cdef int firstA = informative[0].breakA
    cdef int firstB = informative[0].breakB
    cdef int br_a, br_b, seperation
    cdef bint try_cluster = False

    if insert_size != -1 and len(informative) > 1:  # paired end mode
        for v_item in informative:
            br_a = v_item.breakA
            br_b = v_item.breakB
            coords[idx, 0] = br_a
            coords[idx, 1] = br_b
            idx += 1
            if idx > 0 and not try_cluster:
                seperation = np.sqrt((br_a - firstA)**2 + (br_b - firstB)**2)
                if seperation > insert_size:
                    try_cluster = True

    sub_cluster_calls = []
    cdef EventResult_t er
    if try_cluster:
        try:
            Z = linkage(coords, 'single')
        except:
            echo(informative[0].chrA)
            echo(informative[0].breakA)
            echo(coords)
            echo(linkage(coords, 'single'))
            quit()
        clusters = fcluster(Z, insert_size, criterion='distance')
        clusters_d = defaultdict(list)
        for cluster_id, v_item in zip(clusters, informative):
            clusters_d[cluster_id].append(v_item)
        blocks = clusters_d.values()
        for sub_informative in blocks:
            er = make_single_call(sub_informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                    spanning_alignments, extended_tags, 1, generic_insertions)
            sub_cluster_calls.append(er)
    else:
        er = make_single_call(informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                spanning_alignments, extended_tags, 1, generic_insertions)

        sub_cluster_calls.append(er)

    return sub_cluster_calls


cdef single(infile, rds, int insert_size, int insert_stdev, float insert_ppf, int clip_length, int min_support,
                 int to_assemble, int extended_tags):

    # Infer the other breakpoint from a single group
    # Make sure at least one read is worth calling
    # The group may need to be split into multiple calls using the partition_single function
    cdef int min_distance = insert_size + (2*insert_stdev)
    cdef int n_templates = len(set([i.qname for _, i in rds]))

    if n_templates == 1:
        # Filter uninteresting reads
        if not any(not i.flag & 1 or not i.flag & 2 or i.rname != i.rnext or node_info[5] != 2 or
                   (i.flag & 1 and abs(i.tlen) > min_distance)
                   for node_info, i in rds):
            return []

    # Group by template name to check for split reads
    informative = []  # make call from these
    spanning_alignments = []  # make call from these if available (and count attributes)
    generic_insertions = []
    tmp = defaultdict(list)  # group by template name
    small_tlen_outliers = 0  # for paired reads note any smaller than expected TLENs

    for cigar_info, align in rds:
        tmp[align.qname].append((cigar_info, align))

        if align.flag & 1 and abs(align.tlen) < insert_ppf:
            small_tlen_outliers += 1

    cdef AlignmentItem v_item, itm
    cdef int left_clip_a, right_clip_a, left_clip_b, right_clip_b

    for temp_name, alignments in tmp.items():

        l_align = list(alignments)

        if len(l_align) > 1:
            pair = guess_informative_pair(l_align)
            if pair is not None:
                if len(pair) == 2:
                    a, b = pair
                    v_item = make_sv_alignment_item(a, b)
                    classify_d(v_item)
                    informative.append(v_item)

                else:
                    spanning_alignments.append(pair)

        else:  # Single alignment, check spanning
            cigar_info, a = alignments[0]
            cigar_index = cigar_info[5]
            if 0 < cigar_index < len(a.cigartuples) - 1:  # Alignment spans SV
                event_pos = cigar_info[6]
                ci = a.cigartuples[cigar_index]
                spanning_alignments.append((ci[0],
                                            a.rname,
                                            event_pos,
                                            within_read_end_position(event_pos, ci[0], a.cigartuples, cigar_index),
                                            a.cigartuples[cigar_index][1],
                                            a))

            else:
                gi = make_generic_insertion_item(a, insert_size, insert_stdev)
                if gi is not None:
                    generic_insertions.append(gi)

    cdef int min_found_support = 0
    cdef str svtype, jointype
    cdef bint passed
    cdef EventResult_t er
    if len(spanning_alignments) > 0:

        # make call from spanning alignments if possible
        svtype_m = Counter([i[0] for i in spanning_alignments]).most_common()[0][0]
        spanning_alignments = [i for i in spanning_alignments if i[0] == svtype_m]

        posA_arr = [i[2] for i in spanning_alignments]
        posA = int(np.median(posA_arr))
        posA_95 = int(abs(int(np.percentile(posA_arr, [97.5])) - posA))
        posB_arr = [i[3] for i in spanning_alignments]
        posB = int(np.median(posB_arr))
        posB_95 = int(abs(int(np.percentile(posB_arr, [97.5])) - posB))
        chrom = spanning_alignments[0][1]
        svlen = int(np.mean([i[4] for i in spanning_alignments]))
        ab = abs(posB - posA)

        if svlen > 0:
            jitter = ((posA_95 + posB_95) / 2) / svlen
        else:
            jitter = -1

        er = EventResult()
        er.svtype = "DEL" if svtype_m == 2 else "INS"
        er.join_type = "3to5"
        er.chrA = chrom
        er.chrB = chrom
        er.posA = posA
        er.posB = posB
        er.cipos95A = posA_95
        er.cipos95B = posB_95
        er.preciseA = True
        er.preciseB = True
        er.svlen = svlen if svlen >= ab else ab
        er.query_gap = 0
        er.query_overlap = 0
        er.jitter = jitter

        u_reads = [i[5] for i in spanning_alignments]
        v_reads = []
        min_found_support = len(spanning_alignments)

        if len(generic_insertions) > 0:
            min_found_support += len(generic_insertions)
            informative += generic_insertions

        informative_reads = []
        for v_item in informative:
            informative_reads.append(v_item.read_a)
            if v_item.read_b is not None:
                informative_reads.append(v_item.read_b)

        count_attributes2(informative_reads, [], [i[5] for i in spanning_alignments], extended_tags, insert_ppf, generic_insertions, er)

        er.contig = None
        er.contig_left_weight = 0
        er.contig_right_weight = 0
        er.contig2 = None
        er.contig2_left_weight = 0
        er.contig2_right_weight = 0

        as1 = None
        as2 = None
        ref_bases = 0

        if to_assemble: # or len(spanning_alignments) > 0:
            if er.preciseA:
                as1 = assembler.base_assemble(u_reads, er.posA, 500)
                ref_bases += assign_contig_to_break(as1, er, "A", spanning_alignments)
            if er.preciseB:
                as2 = assembler.base_assemble(v_reads, er.posB, 500)
                ref_bases += assign_contig_to_break(as2, er, "B", 0)
            if not as1 and len(generic_insertions) > 0:
                as1 = assembler.base_assemble([item.read_a for item in generic_insertions], er.posA, 500)
                ref_bases += assign_contig_to_break(as1, er, "A", 0)

        er.linked = 0
        er.block_edge = 0
        er.ref_bases = ref_bases
        er.sqc = -1  # not calculated

        return er

    elif len(informative) > 0:
        return partition_single(informative, insert_size, insert_stdev, insert_ppf, spanning_alignments,
                                min_support, extended_tags, to_assemble, [])

    elif len(generic_insertions) > 0:
        # this is a bit confusing, generic insertions are used instead to make call, but bnd count is now 0
        return make_single_call([], insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                spanning_alignments, extended_tags, 0, generic_insertions)


cdef tuple informative_pair(u, v):

    pri_u = None
    pri_v = None
    for i in u:
        ri_flag = i.flag & 64
        if not i.flag & 2304:  # Not pri, supplementary --> is primary
            pri_u = i
        for j in v:
            if j.flag & 64 == ri_flag:  # Same read
                # Same read, primary + supp, or supp + supp
                return i, j

            if not j.flag & 2304:  # Is primary
                pri_v = j

    if pri_u is not None and pri_v is not None:
        return pri_u, pri_v


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

                v.svtype = "INS"
                # Check if gap on query is bigger than gap on reference
                # query_gap = abs(v.b_qstart - v.a_qend) + 1  # as A is first
                ref_gap = abs(v.breakB - v.breakA)

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

    # B is first
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

            # v.svtype = "DEL"
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
            # query_gap = v.b_qstart - v.a_qend  # as A is first
        else:
            query_gap = v.b_qstart - v.a_qend
        if query_gap < 0:
            query_overlap = abs(query_gap)
            query_gap = 0

    v.query_gap = query_gap
    v.query_overlap = query_overlap


cdef void classify_d(AlignmentItem v):

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


cdef tuple break_ops(positions, precise, int limit, float median_pos):

    ops = []
    cdef int i
    for i in positions:
        ops.append((i, 1, 0))
        ops.append((i + limit, -1, 0))
    cdef int v = 1
    if limit < 0:
        v = -1
    for i in precise:
        ops.append((i, 1, 1))
        if v == -1:
            ops.append((i + v + 1, -1, 1))
        else:
            ops.append((i + v, -1, 1))

    ops = sorted(ops, reverse=limit < 0)
    cdef int cum_sum = 0
    cdef int max_sum = 0
    cdef int max_sum_i = 0
    cdef int idx
    for idx, item in enumerate(ops):
        cum_sum += item[1]
        if cum_sum > max_sum:
            max_sum = cum_sum
            max_sum_i = idx

    cdef int break_point, cipos95, is_precise
    break_point, _, is_precise = ops[max_sum_i]
    cipos95 = 0
    if not is_precise:
        if any(i == break_point for i in precise):
            is_precise = 1

    if is_precise and len(precise) > 1:
        cipos95 = 0
    else:
        cipos95 = int(abs(int(np.percentile(positions, [97.5])) - median_pos))

    return break_point, cipos95, is_precise


cdef infer_unmapped_insertion_break_point(int main_A_break, int cipos95A, int preciseA, int main_B_break, int cipos95B,
                                          int preciseB):
    # A or B break pay be in accurate, try and infer svlen and a more accurate break point(s)
    cdef int svlen = 0
    cdef int mid
    if not preciseA and not preciseB:
        # Make break at mid point
        sep = abs(main_B_break - main_A_break)
        svlen = max(sep, max(cipos95A, cipos95B))
        mid = int(sep / 2) + min(main_A_break, main_B_break)
        main_A_break = mid
        main_B_break = mid + 1

    if not preciseA and preciseB:
        if main_A_break < main_B_break:
            main_A_break = main_B_break - 1
        else:
            main_A_break = main_B_break + 1
        svlen = int(cipos95A / 2)

    if not preciseB and preciseA:
        if main_A_break < main_B_break:
            main_B_break = main_A_break + 1
        else:
            main_B_break = main_A_break - 1
        svlen = int(cipos95B / 2)

    return main_A_break, cipos95A, preciseA, main_B_break, cipos95B, preciseB, svlen


cdef void make_call(informative, breakA_precise, breakB_precise, svtype, jointype,
                    int insert_size, int insert_stdev, EventResult_t er):
    # Inspired by mosdepth algorithm +1 for start -1 for end using intervals where break site could occur
    # Use insert size to set a limit on where break site could occur
    cdef int limit = insert_size + insert_stdev
    cdef AlignmentItem i
    # get bulk call
    positionsA = [i.breakA for i in informative]
    positionsB = [i.breakB for i in informative]
    cdef float median_A = np.median(positionsA)
    cdef float median_B = np.median(positionsB)
    cdef int main_A_break = 0
    cdef int main_B_break = 0
    cdef int cipos95A, cipos95B

    if svtype == "DEL":
        if median_A < median_B:
            main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, limit, median_A)
            main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, -limit, median_B)
        else:
            main_B_break, cipos95B, preciseA = break_ops(positionsB, breakB_precise, limit, median_B)
            main_A_break, cipos95A, preciseB = break_ops(positionsA, breakA_precise, -limit, median_A)

    elif svtype == "DUP":
        if median_A < median_B:
            main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, -limit, median_A)
            main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, limit, median_B)
        else:
            main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, limit, median_A)
            main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, -limit, median_B)

    # Generic types
    elif jointype == "3to3":
        main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, limit, median_A)
        main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, limit, median_B)

    elif jointype == "5to5":
        main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, -limit, median_A)
        main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, -limit, median_B)

    elif jointype == "3to5":
        main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, limit, median_A)
        main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, -limit, median_B)

    else:  # 5to3
        main_A_break, cipos95A, preciseA = break_ops(positionsA, breakA_precise, -limit, median_A)
        main_B_break, cipos95B, preciseB = break_ops(positionsB, breakB_precise, limit, median_B)

    svlen_precise = 1
    if informative[0].chrA == informative[0].chrB:
        if svtype == "INS":  # use inferred query gap to call svlen
            svlen_precise = 0
            lens = []
            inferred_lens = []
            for i in informative:
                if i.inferred_sv_len != -1:
                    if i.size_inferred == 1:
                        inferred_lens.append(i.inferred_sv_len)
                    else:
                        lens.append(i.inferred_sv_len)
            if len(lens) > 0:
                svlen = int(np.mean(lens))
                svlen_precise = 1
            elif len(inferred_lens) > 0:
                svlen = int(np.mean(inferred_lens))
            else:
                main_A_break, cipos95A, preciseA, main_B_break, cipos95B, preciseB, svlen = \
                    infer_unmapped_insertion_break_point(main_A_break, cipos95A, preciseA, main_B_break, cipos95B, preciseB)
        elif svtype == "DEL":

            main_svlen = abs(main_B_break - main_A_break)
            svlen = main_svlen
            if not preciseA or not preciseB:

                lens = []
                inferred_lens = []
                for i in informative:
                    if i.inferred_sv_len != -1:
                        if i.size_inferred == 1:
                            inferred_lens.append(i.inferred_sv_len)
                        else:
                            lens.append(i.inferred_sv_len)
                if len(lens) > 0:
                    svlen = int(np.mean(lens))
                    if main_svlen > 0 and (svlen / main_svlen) > 0.7:
                        svlen_precise = 1
                    else:
                        svlen = main_svlen
                else:
                    if len(inferred_lens) > 0:
                        svlen = int(np.mean(inferred_lens))
                    else:
                        svlen = main_svlen

                    if main_svlen > 0 and (svlen / main_svlen) <= 0.7:
                        svlen = main_svlen

        else:
            svlen = abs(main_B_break - main_A_break)
    else:
        svlen = 0

    q_overlaps = int(np.mean([i.query_overlap for i in informative]))
    if q_overlaps != q_overlaps:
        q_overlaps = 0

    if svlen > 0:
        jitter = ((cipos95A + cipos95B) / 2) / svlen
    else:
        jitter = -1

    er.svtype = svtype
    er.join_type = jointype
    er.chrA = informative[0].chrA
    er.chrB = informative[0].chrB

    er.cipos95A = cipos95A
    er.cipos95B = cipos95B
    er.posA = main_A_break
    er.posB = main_B_break
    er.preciseA = preciseA
    er.preciseB = preciseB
    er.svlen_precise = svlen_precise
    er.svlen = svlen
    er.query_overlap = q_overlaps
    er.jitter = jitter


cdef tuple mask_soft_clips(int aflag, int bflag, a_ct, b_ct):

    # Find out which soft clip pairs are compatible with the chosen read pair
    cdef int left_clipA = 1 if (a_ct[0][0] == 4 or a_ct[0][0] == 5) else 0
    cdef int right_clipA = 1 if (a_ct[-1][0] == 4 or a_ct[-1][0] == 5) else 0
    cdef int left_clipB = 1 if (b_ct[0][0] == 4 or b_ct[0][0] == 5) else 0
    cdef int right_clipB = 1 if (b_ct[-1][0] == 4 or b_ct[-1][0] == 5) else 0

    cdef int a_template_start = 0
    cdef int b_template_start = 0

    if aflag & 64 == bflag & 64:  # Same read

        if (left_clipB and right_clipB) or (left_clipA and right_clipA):  # One read has more than one soft-clip

            if aflag & 16:  # A is reverse, convert to forward
                if a_ct[-1][0] != 0:
                    a_template_start = a_ct[-1][1]
            else:
                if a_ct[0][0] != 0:
                    a_template_start = a_ct[0][1]
            if bflag & 16:  # B is reverse
                if b_ct[-1][0] != 0:
                    b_template_start = b_ct[-1][1]
            else:
                if b_ct[0][0] != 0:
                    b_template_start = b_ct[0][1]

            if left_clipB and right_clipB:  # Choose one soft-clip for B

                if b_template_start < a_template_start:
                    if not bflag & 16:  # B on forward strand
                        left_clipB = 0
                    else:
                        right_clipB = 0
                else:
                    if not bflag & 16:  # B on reverse strand
                        right_clipB = 0
                    else:
                        left_clipB = 0

            if left_clipA and right_clipA:  # Choose one soft-clip for A

                if a_template_start < b_template_start:
                    if not aflag & 16:  # A on forward strand
                        left_clipA = 0
                    else:
                        right_clipA = 0
                else:
                    if not aflag & 16:  # A on reverse strand
                        right_clipA = 0
                    else:
                        left_clipA = 0

    else:  # Different reads choose longest if more than one soft-clip
        if left_clipB and right_clipB:
            if b_ct[0][1] > b_ct[-1][1]:
                right_clipB = 0
            else:
                left_clipB = 0

        if left_clipA and right_clipA:
            if a_ct[0][1] > a_ct[-1][1]:
                right_clipA = 0
            else:
                left_clipA = 0

    return left_clipA, right_clipA, left_clipB, right_clipB


cdef query_start_end_from_cigartuples(r):
    # Infer the position on the query sequence of the alignment using cigar string
    cdef int end = 0
    cdef int start = 0
    cdef bint i = 0
    cdef int opp, slen
    for opp, slen in r.cigartuples:
        if opp == 0:
            end += slen
        elif opp == 4 or opp == 5:
            if not i:
                start += slen
                end += slen
                i = 1
            else:
                break
        elif opp == 1:
            end += slen
        i = 1
    return start, end


cdef start_end_query_pair(r1, r2):
    cdef int query_length, s1, e1, s2, e2, start_temp
    # r1 and r2 might be on same read, if this is case infer the query position on the read
    s1, e1 = query_start_end_from_cigartuples(r1)
    s2, e2 = query_start_end_from_cigartuples(r2)
    if r1.flag & 64 == r2.flag & 64:  # same read
        if r2.flag & 16 != r1.flag & 16:  # different strand, count from end
            query_length = r1.infer_read_length()  # Note, this also counts hard-clips
            start_temp = query_length - e2
            e2 = start_temp + e2 - s2
            s2 = start_temp
    return s1, e1, s2, e2


def sort_by_length(x):
    return len(x)


cdef void assemble_partitioned_reads(EventResult_t er, u_reads, v_reads, int block_edge, int assemble):

    as1 = None
    as2 = None
    if assemble:
        if er.preciseA:
            as1 = assembler.base_assemble(u_reads, er.posA, 500)
            if as1:
                if er.spanning == 0 and not (as1['left_clips'] or as1['right_clips']):
                    as1 = None
        if (er.spanning == 0 or as1 is None) and er.preciseB:
            as2 = assembler.base_assemble(v_reads, er.posB, 500)
            if as2 :
                if not (as2['left_clips'] or as2['right_clips']):
                    as2 = None

    er.linked = 0
    er.block_edge = block_edge
    er.contig = None
    er.contig_left_weight = 0
    er.contig_ref_start = -1
    er.contig_ref_end = -1
    er.contig_right_weight = 0
    er.contig2 = None
    er.contig2_left_weight = 0
    er.contig2_right_weight = 0
    er.contig2_ref_start = -1
    er.contig2_ref_end = -1
    rbases = 0
    if as1 is not None and "contig" in as1:
        er.contig = as1["contig"]
        rbases += as1["ref_bases"]
        er.contig_ref_start = as1["ref_start"]
        er.contig_ref_end = as1["ref_end"]
        er.contig_left_weight = as1["left_weight"]
        er.contig_right_weight = as1["right_weight"]
    if as2 is not None and "contig" in as2:
        er.contig2 = as2["contig"]
        rbases += as2["ref_bases"]
        er.contig2_ref_start = as2["ref_start"]
        er.contig2_ref_end = as2["ref_end"]
        er.contig2_left_weight = as2["left_weight"]
        er.contig2_right_weight = as2["right_weight"]
    er.ref_bases = rbases


cdef call_from_reads(u_reads, v_reads, int insert_size, int insert_stdev, float insert_ppf, int min_support, int block_edge, int assemble, int extended_tags):

    grp_u = defaultdict(list)
    grp_v = defaultdict(list)

    for r in u_reads:
        grp_u[r.qname].append(r)

    for r in v_reads:
        grp_v[r.qname].append(r)

    informative = []

    cdef AlignmentItem v_item, i
    cdef int left_clip_a, right_clip_a, left_clip_b, right_clip_b

    for qname in [k for k in grp_u if k in grp_v]:  # Qname found on both sides

        u = grp_u[qname]
        v = grp_v[qname]

        if len(u) == 1 and len(v) == 1:
            pair = (u[0], v[0])
        else:
            pair = informative_pair(u, v)

        if pair:

            a, b = pair
            a_ct = a.cigartuples
            b_ct = b.cigartuples

            a_qstart, a_qend, b_qstart, b_qend = start_end_query_pair(a, b)

            # Soft-clips for the chosen pair, plus template start of alignment
            left_clip_a, right_clip_a, left_clip_b, right_clip_b = mask_soft_clips(a.flag, b.flag, a_ct, b_ct)

            a_chrom, b_chrom = a.rname, b.rname
            a_start, a_end = a.pos, a.reference_end
            b_start, b_end = b.pos, b.reference_end
            read_overlaps_mate = same_read_overlaps_mate(a_chrom, b_chrom, a_start, a_end, b_start, b_end, a, b)

            v_item = AlignmentItem(a.rname,
                                   b.rname,
                                   int(not a.flag & 2304),  # Is primary
                                   int(not b.flag & 2304),
                                   1 if a.flag & 64 else 2,
                                   1 if b.flag & 64 else 2,
                                   a.pos,
                                   a.reference_end,
                                   b.pos,
                                   b.reference_end,
                                   3 if a.flag & 16 == 0 else 5,
                                   3 if b.flag & 16 == 0 else 5,
                                   left_clip_a,
                                   right_clip_a,
                                   left_clip_b,
                                   right_clip_b,
                                   a_qstart, a_qend, b_qstart, b_qend,
                                   read_overlaps_mate,
                                   a, b
                                   )

            if v_item.left_clipA and v_item.right_clipA:
                if a_ct[0][1] >= a_ct[-1][1]:
                    v_item.right_clipA = 0
                else:
                    v_item.left_clipA = 0

            if v_item.left_clipB and v_item.right_clipB:
                if b_ct[0][1] >= b_ct[-1][1]:
                    v_item.right_clipB = 0
                else:
                    v_item.left_clipB = 0

            classify_d(v_item)
            informative.append(v_item)

    if not informative:
        return []

    svtypes_counts = [[], [], [], []]
    for i in informative:
        if i.svtype == "DEL":
            svtypes_counts[0].append(i)
        elif i.svtype == "DUP":
            svtypes_counts[1].append(i)
        elif i.svtype == "INV":
            svtypes_counts[2].append(i)
        elif i.svtype == "TRA":
            svtypes_counts[3].append(i)

    svtypes_res = sorted(svtypes_counts, key=sort_by_length, reverse=True)

    results = []
    cdef EventResult_t er
    for sub_informative in svtypes_res:
        if len(sub_informative) >= min_support:
            svtype = sub_informative[0].svtype
            if svtype == "INV" or svtype == "TRA":
                jointype = Counter([i.join_type for i in sub_informative]).most_common(1)[0][0]
            else:
                jointype = sub_informative[0].join_type

            precise_a = []
            precise_b = []
            u_reads = []
            v_reads = []

            for v_item in sub_informative:
                if v_item.breakA_precise:
                    precise_a.append(v_item.breakA)
                if v_item.breakB_precise:
                    precise_b.append(v_item.breakB)
                u_reads.append(v_item.read_a)
                v_reads.append(v_item.read_b)

            er = EventResult()

            make_call(sub_informative, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev, er)

            count_attributes2(u_reads, v_reads, [], extended_tags, insert_ppf, [], er)
            if er.su < min_support:
                continue

            assemble_partitioned_reads(er, u_reads, v_reads, block_edge, assemble)
            results.append(er)

    return results


cdef one_edge(infile, u_reads_info, v_reads_info, int clip_length, int insert_size, int insert_stdev, float insert_ppf,
                   int min_support, int block_edge, int assemble, int extended_tags):

    spanning_alignments = []
    u_reads = []
    v_reads = []
    cdef int cigar_index, event_pos

    for cigar_info, a in u_reads_info:
        if not a.cigartuples:
            continue
        u_reads.append(a)
        cigar_index = cigar_info[5]
        if 0 < cigar_index < len(a.cigartuples) - 1:  # Alignment spans SV
            event_pos = cigar_info[6]
            ci = a.cigartuples[cigar_index]
            spanning_alignments.append((ci[0],
                                        a.rname,
                                        event_pos,
                                        event_pos + 1 if ci[0] == 1 else event_pos + a.cigartuples[cigar_index][1] + 1,
                                        event_pos + a.cigartuples[cigar_index][1] + 1,
                                        a))

    for cigar_info, a in v_reads_info:
        if not a.cigartuples:
            continue
        v_reads.append(a)
        cigar_index = cigar_info[5]
        if 0 < cigar_index < len(a.cigartuples) - 1:  # Alignment spans SV
            event_pos = cigar_info[6]
            ci = a.cigartuples[cigar_index]
            spanning_alignments.append((ci[0],
                                        a.rname,
                                        event_pos,
                                        event_pos + 1 if ci[0] == 1 else event_pos + a.cigartuples[cigar_index][1] + 1,
                                        event_pos + a.cigartuples[cigar_index][1],
                                        a))

    # info = {}
    cdef str svtype, jointype

    if len(spanning_alignments) > 0:
        svtype_m = Counter([i[0] for i in spanning_alignments]).most_common()[0][0]
        spanning_alignments = [i for i in spanning_alignments if i[0] == svtype_m]

    # make call from spanning alignments if possible
    cdef EventResult_t er
    if len(spanning_alignments) > 0:
        posA_arr = [i[2] for i in spanning_alignments]
        posA = int(np.median(posA_arr))
        posA_95 = int(abs(int(np.percentile(posA_arr, [97.5])) - posA))
        posB_arr = [i[3] for i in spanning_alignments]
        posB = int(np.median(posB_arr))
        posB_95 = int(abs(int(np.percentile(posB_arr, [97.5])) - posB))
        chrom = spanning_alignments[0][1]
        svlen = int(np.mean([i[4] for i in spanning_alignments]))

        if svlen > 0:
            jitter = ((posA_95 + posB_95) / 2) / svlen
        else:
            jitter = -1

        er = EventResult()
        er.svtype = "DEL" if svtype_m == 2 else "INS"
        er.join_type = "3to5"
        er.chrA = chrom
        er.chrB = chrom
        er.posA = posA
        er.posB = posB
        er.cipos95A = posA_95
        er.cipos95B = posB_95,
        er.jitter = jitter,
        er.preciseA = True
        er.preciseB = True
        er.svlen = svlen
        er.query_overlap = 0

        u_reads = [i[5] for i in spanning_alignments]
        v_reads = []
        count_attributes2(u_reads, v_reads, [i[5] for i in spanning_alignments], extended_tags, insert_ppf, [], er)
        # if not attrs or attrs["su"] < min_support:
        if er.su < min_support:
            return []

        assemble_partitioned_reads(er, u_reads, v_reads, block_edge, assemble)

        return [er]

    else:
        results = call_from_reads(u_reads, v_reads, insert_size, insert_stdev, insert_ppf, min_support, block_edge, assemble, extended_tags)

        return results


cdef list get_reads(infile, nodes_info, buffered_reads, n2n, bint add_to_buffer):

    cdef int j, int_node
    cdef uint64_t p
    cdef uint64_t v
    cdef AlignedSegment a
    aligns = []
    for int_node in nodes_info:

        n = n2n[int_node]
        if int_node in buffered_reads:
            aligns.append((n, buffered_reads[int_node]))
            continue

        p = n[4]
        node = (n[0], n[1], n[2], n[3], p)  # drop cigar index and event pos
        infile.seek(p)
        try:
            a = next(infile)
        except StopIteration:
            return aligns
        v = xxhasher(bam_get_qname(a._delegate), len(a.qname), 42)

        if (v, a.flag, a.pos, a.rname, p) == node:
            aligns.append((n, a))
            if add_to_buffer:
                buffered_reads[int_node] = a  # Add to buffer, then block nodes with multi-edges dont need collecting twice
            continue
        else:  # Try next few reads, not sure why this occurs
            steps = 0
            while steps < 5:
                try:
                    a = next(infile)
                except StopIteration:
                    return aligns
                steps += 1
                v = xxhasher(bam_get_qname(a._delegate), len(a.qname), 42)
                if (v, a.flag, a.pos, a.rname, p) == node:
                    aligns.append((n, a))
                    if add_to_buffer:
                        buffered_reads[int_node] = a
                    break

    return aligns


cdef list multi(data, bam, int insert_size, int insert_stdev, float insert_ppf, int clip_length, int min_support, int lower_bound_support,
                 int assemble_contigs, int extended_tags):

    # Sometimes partitions are not linked, happens when there is not much support between partitions
    # Then need to decide whether to call from a single partition
    n2n = data["n2n"]
    seen = set(data["parts"].keys())

    out_counts = defaultdict(int)  # The number of 'outward' links to other clusters
    cdef int buffered_reads = 0  # Only buffer reads for smallish components, otherwise memory issues
    cdef bint add_to_buffer = 1
    events = []
    for (u, v), d in data["s_between"].items():

        rd_u = get_reads(bam, d[0], data["reads"], data["n2n"], add_to_buffer)   # [(Nodeinfo, alignment)..]
        rd_v = get_reads(bam, d[1], data["reads"], data["n2n"], add_to_buffer)

        total_reads = len(rd_u) + len(rd_v)

        buffered_reads += total_reads
        if add_to_buffer and buffered_reads > 50000:
            add_to_buffer = 0

        out_counts[u] += len(rd_u)
        out_counts[v] += len(rd_v)
        if len(rd_u) == 0 or len(rd_v) == 0:
            continue

        if u in seen:
            seen.remove(u)
        if v in seen:
            seen.remove(v)

        events += one_edge(bam, rd_u, rd_v, clip_length, insert_size, insert_stdev, insert_ppf, min_support, 1, assemble_contigs,
                               extended_tags)

    # Process any unconnected blocks
    if seen:
        for part_key in seen:
            d = data["parts"][part_key]
            if len(d) >= lower_bound_support:
                # Call single block, only collect local reads to the block
                rds = get_reads(bam, d, data["reads"], data["n2n"], 0)
                if len(rds) < lower_bound_support:
                    continue

                res = single(bam, rds, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs, extended_tags)
                if res:
                    if isinstance(res, EventResult):
                        events.append(res)
                    else:
                        events += res

    # Check for events within clustered nodes
    for k, d in data["s_within"].items():
        o_count = out_counts[k]
        i_counts = len(d)

        if o_count > 0 and i_counts > (2*min_support) and i_counts > o_count:
            rds = get_reads(bam, d, data["reads"], data["n2n"], 0)
            if len(rds) < lower_bound_support:
                continue

            res = single(bam, rds, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs, extended_tags)

            if res:
                if isinstance(res, EventResult):
                    events.append(res)
                else:
                    events += res

    return events


cpdef list call_from_block_model(bam, data, clip_length, insert_size, insert_stdev, insert_ppf, min_support, lower_bound_support, extended_tags,
                                 assemble_contigs):

    n_parts = len(data["parts"])
    events = []

    if n_parts >= 1:
        # Processed single edges and break apart connected
        events += multi(data, bam, insert_size, insert_stdev, insert_ppf, clip_length, min_support, lower_bound_support,
                        assemble_contigs, extended_tags)

    elif n_parts == 0:
        # Possible single read only
        rds = get_reads(bam, data["n2n"].keys(), data["reads"], data["n2n"], 0)

        if len(rds) < lower_bound_support:
            return []

        e = single(bam, rds, insert_size, insert_stdev, insert_ppf, clip_length, min_support,
                   assemble_contigs, extended_tags)
        if e:
            if isinstance(e, list):
                events += e
            else:
                events.append(e)

    return events
