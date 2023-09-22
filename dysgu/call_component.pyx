#cython: language_level=3, boundscheck=True, c_string_type=unicode, c_string_encoding=utf8, infer_types=True
from __future__ import absolute_import
from collections import Counter, defaultdict
import logging
import numpy as np
cimport numpy as np
from numpy.random import normal
import itertools

from dysgu import assembler
from dysgu.map_set_utils import echo
from dysgu.map_set_utils cimport hash as xxhasher
from dysgu.map_set_utils cimport is_overlapping, clip_sizes_hard, EventResult, clip_sizes, min_fractional_overlapping
from dysgu.sv_category cimport AlignmentItem, classify_d
from dysgu.extra_metrics cimport soft_clip_qual_corr
from dysgu.extra_metrics import filter_poorly_aligned_ends, gap_size_upper_bound
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_qname, bam_get_cigar
from libc.stdint cimport uint32_t
from libc.stdint cimport uint64_t
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

np.random.seed(1)

ctypedef EventResult EventResult_t

ctypedef enum ReadEnum_t:
    DISCORDANT = 0
    SPLIT = 1
    DELETION = 2
    INSERTION = 3
    BREAKEND = 4


cdef n_aligned_bases(AlignedSegment aln):
    cdef int opp, l, aligned, large_gaps, n_small_gaps, i
    cdef uint32_t cigar_value
    cdef uint32_t cigar_l = aln._delegate.core.n_cigar
    cdef uint32_t *cigar_p = bam_get_cigar(aln._delegate)
    aligned, large_gaps, n_small_gaps = 0, 0, 0
    for i in range(cigar_l):
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & 15
        l = <int> cigar_value >> 4
        if opp == 0:
            aligned += l
        elif opp == 1 or opp == 2:
            if l >= 30:
                large_gaps += l
            else:
                n_small_gaps += 1
    return float(aligned), float(large_gaps), float(n_small_gaps)


cdef base_quals_aligned_clipped(AlignedSegment a):
    cdef int left_clip, right_clip
    cdef float aligned_base_quals = 0
    cdef float aligned_bases = 0
    cdef float clipped_base_quals = 0
    left_clip, right_clip = clip_sizes(a)
    clipped_bases = left_clip + right_clip
    cdef const unsigned char[:] quals = a.query_qualities
    cdef int i
    for i in range(left_clip):
        clipped_base_quals += quals[i]
    for i in range(left_clip, len(quals) - right_clip):
        aligned_base_quals += quals[i]
        aligned_bases += 1
    for i in range(len(quals) - right_clip, len(quals)):
        clipped_base_quals += quals[i]
    return aligned_base_quals, aligned_bases, clipped_base_quals, clipped_bases


cdef count_attributes2(reads1, reads2, spanning, float insert_ppf, generic_ins,
                       EventResult_t er, bint paired_end_reads):
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
    cdef float aligned_base_quals = 0
    cdef float aligned_bases = 0
    cdef float clipped_base_quals = 0
    cdef float clipped_bases = 0
    cdef int abq, ab, cbq, cb
    paired_end = set([])
    seen = set([])
    er.spanning = len(spanning)
    er.bnd = len(generic_ins)
    cdef int flag, index, this_as
    cdef float a_bases, large_gaps, n_small_gaps
    cdef AlignedSegment a
    for index, a in enumerate(itertools.chain(reads1, reads2, [i.read_a for i in generic_ins])):
        flag = a.flag
        if flag & 2:
            er.NP += 1
        if a.flag & 1 and a.tlen and abs(a.tlen) < insert_ppf:
            er.n_small_tlen += 1
        if paired_end_reads and paired_end and flag & 8:
            er.n_unmapped_mates += 1
        left_clip, right_clip = clip_sizes_hard(a)
        if left_clip > 0 and right_clip > 0:
            er.double_clips += 1
        has_sa = a.has_tag("SA")
        if has_sa:
            n_sa += a.get_tag("SA").count(";")
        if a.has_tag("XA"):
            n_xa += a.get_tag("XA").count(";")
        a_bases, large_gaps, n_small_gaps = n_aligned_bases(a)
        if a_bases > 0:
            n_gaps += n_small_gaps / a_bases
        if flag & 2304:  # Supplementary (and not primary if -M if flagged using bwa)
            er.supp += 1
            MAPQsupp += a.mapq
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
            if paired_end_reads:
                if index >= len(reads2) and a.qname in paired_end:
                    er.pe += 1
                else:
                    paired_end.add(a.qname)
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
        if a.flag & 1:  # paired read
            abq, ab, cbq, cb = base_quals_aligned_clipped(a)
            aligned_base_quals += abq
            aligned_bases += ab
            clipped_base_quals += cbq
            clipped_bases += cb

    for a in spanning:
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
        a_bases, large_gaps, n_small_gaps = n_aligned_bases(a)
        if a_bases > 0:
            n_gaps += n_small_gaps / a_bases
        if flag & 2304:  # Supplementary
            MAPQsupp += a.mapq
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
            if paired_end_reads:
                if a.qname in paired_end:  # If two primary reads from same pair
                    er.pe += 2
                else:
                    paired_end.add(a.qname)
            if a.has_tag("NM"):
                if a_bases:
                    nm = float(a.get_tag("NM"))
                    NMpri += nm / a_bases
                    NMbase += (nm - large_gaps) / a_bases
        if flag & 16:
            er.minus += 1
        else:
            er.plus += 1
        if a.flag & 1:
            abq, ab, cbq, cb = base_quals_aligned_clipped(a)
            aligned_base_quals += abq
            aligned_bases += ab
            clipped_base_quals += cbq
            clipped_bases += cb

    cdef int tot = er.supp + total_pri
    er.NMpri = (NMpri / total_pri) * 100 if total_pri > 0 else 0
    er.NMsupp = (NMsupp / er.supp) * 100 if er.supp > 0 else 0
    er.NMbase = (NMbase / tot) * 100 if tot > 0 else 0
    er.n_gaps = (n_gaps / tot) * 100 if tot > 0 else 0
    er.MAPQpri = MAPQpri / total_pri if total_pri > 0 else 0
    er.MAPQsupp = MAPQsupp / er.supp if er.supp > 0 else 0
    er.n_sa = n_sa / tot if tot > 0 else 0
    er.n_xa = n_xa / tot if tot > 0 else 0
    er.maxASsupp = maxASsupp
    if paired_end_reads and len(reads2) == 0:
        er.pe = len(reads1)
    er.su = er.pe + er.supp + (2*er.spanning) + er.bnd
    if clipped_bases > 0 and aligned_bases > 0 and clipped_base_quals > 0:
        er.clip_qual_ratio = (aligned_base_quals / aligned_bases) / (clipped_base_quals / clipped_bases)
    else:
        er.clip_qual_ratio = 0


cdef int within_read_end_position(event_pos, svtype, cigartuples, cigar_index):
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
        if a.flag & 1 and a_cigar_info.cigar_index == -1 and b_cigar_info.cigar_index == -1:
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
        if 0 < a_cigar_info.cigar_index < len(a.cigartuples) - 1:
            cigar_index = a_cigar_info.cigar_index
            event_pos = a_cigar_info.event_pos
            ci = a.cigartuples[cigar_index]
            return (ci[0],
                    a.rname,
                    event_pos,
                    event_pos + 1 if ci[0] == 1 else event_pos + a.cigartuples[cigar_index][1] + 1,
                    a.cigartuples[cigar_index][1],
                    a,
                    cigar_index)
        elif 0 < b_cigar_info.cigar_index < len(b.cigartuples) - 1:
            cigar_index = b_cigar_info.cigar_index
            event_pos = b_cigar_info.event_pos
            ci = b.cigartuples[cigar_index]
            return (ci[0],
                    b.rname,
                    event_pos,
                    event_pos + 1 if ci[0] == 1 else event_pos + b.cigartuples[cigar_index][1] + 1,
                    b.cigartuples[cigar_index][1],
                    b,
                    cigar_index)
        a_pos = a_cigar_info.event_pos  # Position may have been inferred from SA tag, use this if available
        b_pos = b_cigar_info.event_pos
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
    a_qstart, a_qend, b_qstart, b_qend, a_len, b_len = start_end_query_pair(a, b)
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
                           a_qstart, a_qend, b_qstart, b_qend, a_len, b_len,
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
    aln_span = aln.reference_end - aln.pos
    v_item.size_inferred = 1
    if insert_std > 0:
        rand_insert_pos = abs(insert_size - aln_span + int(normal(0, insert_std)))
    else:  # single read mode
        v_item.svtype = "BND"
        clip_s = max(clip_sizes(aln))
        rand_insert_pos = 100 if not clip_s else clip_s
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
                er.contig2_lc = asmb["left_clips"]
                er.contig2_rc = asmb["right_clips"]
            else:
                er.contig = asmb["contig"]
                er.contig_ref_start = asmb["ref_start"]
                er.contig_ref_end = asmb["ref_end"]
                er.contig_left_weight = asmb["left_weight"]
                er.contig_right_weight = asmb["right_weight"]
                er.contig_lc = asmb["left_clips"]
                er.contig_rc = asmb["right_clips"]
        else:
            # assign contigs to current side
            if side == "A":
                er.contig = asmb["contig"]
                er.contig_ref_start = asmb["ref_start"]
                er.contig_ref_end = asmb["ref_end"]
                er.contig_left_weight = asmb["left_weight"]
                er.contig_right_weight = asmb["right_weight"]
                er.contig_lc = asmb["left_clips"]
                er.contig_rc = asmb["right_clips"]
            else:
                er.contig2 = asmb["contig"]
                er.contig2_ref_start = asmb["ref_start"]
                er.contig2_ref_end = asmb["ref_end"]
                er.contig2_left_weight = asmb["left_weight"]
                er.contig2_right_weight = asmb["right_weight"]
                er.contig2_lc = asmb["left_clips"]
                er.contig2_rc = asmb["right_clips"]
        ref_bases += asmb["ref_bases"]
    return ref_bases


cdef make_single_call(sub_informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble, spanning_alignments,
                      svlen_precise, generic_ins, site, bint paired_end):
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
        if v_item.read_a is not None:
            u_reads.append(v_item.read_a)
        if v_item.read_b is not None:
            v_reads.append(v_item.read_b)
    if not u_reads and not v_reads:
        return

    call_informative = Counter([(itm.svtype, itm.join_type) for itm in si]).most_common()
    svtype, jointype = call_informative[0][0]
    make_call(si, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev, er)
    if len(sub_informative) > 0:
        count_attributes2(u_reads, v_reads, [], insert_ppf, generic_ins, er, paired_end)
    else:
        count_attributes2([], [], [], insert_ppf, generic_ins, er, paired_end)

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
            if as1 and (er.svtype != "TRA" or (as1['contig'] and (as1['contig'][0].islower() or as1['contig'][-1].islower()))):
                ref_bases += assign_contig_to_break(as1, er, "A", spanning_alignments)
        if er.preciseB:
            as2 = assembler.base_assemble(v_reads, er.posB, 500)
            if as2 and (er.svtype != "TRA" or (as2['contig'] and (as2['contig'][0].islower() or as2['contig'][-1].islower()))):
                ref_bases += assign_contig_to_break(as2, er, "B", 0)
    er.linked = 0
    er.block_edge = 0
    er.ref_bases = ref_bases
    er.svlen_precise = svlen_precise  # if 0 then soft-clip will be remapped
    if svlen_precise == 0:
        corr_score = soft_clip_qual_corr(u_reads + v_reads)
        er.sqc = corr_score
    else:
        er.sqc = -1
    if site:
        er.site_info = site
    return er


def assign_sites_to_clusters(sites_info, clusters, informative, coords, cluster_count):
    if cluster_count == 1 and len(sites_info) == 1:
        return {1: informative},  {1: sites_info[0]}
    clusters_d = defaultdict(list)
    for cluster_id, v_item in zip(clusters, informative):
        clusters_d[cluster_id].append(v_item)
    # reads have cluster ids, assign any sites to these ids
    # get the mean coords of each cluster
    c_count = np.zeros(cluster_count)
    x_sum = np.zeros_like(c_count)
    y_sum = np.zeros_like(c_count)
    for i in range(len(coords)):
        c_id = clusters[i]
        x_sum[c_id - 1] += coords[i, 0]  # c_id are >= 1
        y_sum[c_id - 1] += coords[i, 1]
        c_count[c_id - 1] += 1
    m_x = x_sum / c_count
    m_y = x_sum / c_count
    c_xy = {}
    for cluster_id in clusters_d.keys():
        c_xy[cluster_id] = (m_x[cluster_id - 1], m_y[cluster_id - 1])

    # assign sites to clusters in greedy fashion
    sites_to_clusters = {}
    un_partitoned = sites_info
    c = cluster_count
    while len(un_partitoned) > 0 and c > 0 and len(c_xy) > 0:
        s = un_partitoned.pop()
        best_id = -1
        best_dist = 1e9
        for cluster_id, (x, y) in c_xy.items():
            sep = np.sqrt((s.start - x)**2 + (s.end - y)**2)
            if sep < best_dist:
                best_id = cluster_id
                best_dist = sep
        if best_id == -1:
            raise ValueError("best_id", best_id, coords, sites_info, cluster_id)
        sites_to_clusters[best_id] = s
        del c_xy[best_id]
        c -= 1
    return clusters_d, sites_to_clusters


cdef partition_single(informative, insert_size, insert_stdev, insert_ppf, spanning_alignments,
                      min_support, to_assemble, generic_insertions, sites_info, bint paired_end):
    # spanning alignments is empty
    cdef AlignmentItem v_item
    cdef int idx = 0
    cdef np.ndarray[double, ndim=2] coords = np.zeros((len(informative) + len(sites_info), 2))
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
                sep = np.sqrt((br_a - firstA)**2 + (br_b - firstB)**2)
                if sep > insert_size:
                    try_cluster = True
        # sites added here to influence clustering
        if sites_info and len(informative) > 1:
            try_cluster = True
            for ii in sites_info:
                coords[idx, 0] = ii.start
                coords[idx, 1] = ii.end
                idx += 1
    sub_cluster_calls = []
    cdef EventResult_t er
    if try_cluster:
        try:
            Z = linkage(coords, 'single')
        except:
            logging.warning("Linkage clustering failed, at chromosome tid {}, position {}. n-reads in cluster {}".format(
                informative[0].chrA, informative[0].breakA, len(coords)))
            return []
        try:
            clusters = fcluster(Z, insert_size, criterion='distance')
        except:
            logging.warning("fcluster function failed, at chromosome tid {}, position {}. n-reads in cluster {}".format(
                informative[0].chrA, informative[0].breakA, len(coords)))
            return []
        # cluster ids start are >=1, so bincount 0 is always 0
        cluster_count = len(np.bincount(clusters)) - 1
        if sites_info:
            # assign sites to exactly one cluster
            clusters_d, sites_2_clusters = assign_sites_to_clusters(sites_info, clusters, informative, coords, cluster_count)
            for cid, sub_informative in clusters_d.items():
                if cid in sites_2_clusters:
                    st = sites_2_clusters[cid]
                else:
                    st = None
                er = make_single_call(sub_informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                      spanning_alignments, 1, generic_insertions, st, paired_end)
                sub_cluster_calls.append(er)
        else:
            if cluster_count == 1:
                er = make_single_call(informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                      spanning_alignments, 1, generic_insertions, None, paired_end)
                sub_cluster_calls.append(er)
            else:
                clusters_d = defaultdict(list)
                for cluster_id, v_item in zip(clusters, informative):
                    clusters_d[cluster_id].append(v_item)
                for sub_informative in clusters_d.values():
                    er = make_single_call(sub_informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                          spanning_alignments, 1, generic_insertions, None, paired_end)
                    sub_cluster_calls.append(er)
    else:
        if sites_info:
            st = sites_info[0]
        else:
            st = None
        er = make_single_call(informative, insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                              spanning_alignments, 1, generic_insertions, st, paired_end)
        sub_cluster_calls.append(er)
    return sub_cluster_calls


cdef single(rds, int insert_size, int insert_stdev, float insert_ppf, int clip_length, int min_support,
            int to_assemble, sites_info, bint paired_end, int length_extend, float divergence):

    # Infer the other breakpoint from a single group
    # Make sure at least one read is worth calling
    # The group may need to be split into multiple calls using the partition_single function
    cdef int min_distance = insert_size + (2*insert_stdev)
    cdef int n_templates = len(set([i.qname for _, i in rds]))
    if n_templates == 1:
        if not any(not i.flag & 1 or not i.flag & 2 or i.rname != i.rnext or node_info.cigar_index != 2 or
                   (i.flag & 1 and abs(i.tlen) > min_distance)
                   for node_info, i in rds):
            return []
    # Group by template name to check for split reads
    informative = []  # make call from these
    spanning_alignments = []  # make call from these, if available (and count attributes)
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
            cigar_index = cigar_info.cigar_index

            if 0 < cigar_index < len(a.cigartuples) - 1:  # Alignment spans SV
                event_pos = cigar_info.event_pos
                ci = a.cigartuples[cigar_index]
                spanning_alignments.append((ci[0],
                                            a.rname,
                                            event_pos,
                                            within_read_end_position(event_pos, ci[0], a.cigartuples, cigar_index),
                                            a.cigartuples[cigar_index][1],
                                            a,
                                            cigar_index))
            else:
                gi = make_generic_insertion_item(a, insert_size, insert_stdev)
                if gi is not None:
                    generic_insertions.append(gi)

    cdef int min_found_support = 0
    cdef str svtype, jointype
    cdef bint passed
    cdef EventResult_t er
    if len(spanning_alignments) > 0:

        if not paired_end:
            spanning_alignments, rate_poor_ends = filter_poorly_aligned_ends(spanning_alignments, divergence)
            if not spanning_alignments or rate_poor_ends > 0.7:
                return []

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
        # choose representative alignment to use
        best_index = 0
        best_dist = 1e9
        for index in range(len(spanning_alignments)):
            dist = abs(spanning_alignments[index][2] - posA) + abs(spanning_alignments[index][3] - posB)
            if dist < best_dist:
                best_index = index
                best_dist = dist
                if dist == 0:
                    break

        best_align = spanning_alignments[best_index][5]

        er = EventResult()

        if not paired_end:
            size_pos_bounded = [gap_size_upper_bound(sp[5], sp[6], sp[2], sp[3], length_extend, divergence) for sp in spanning_alignments]
            svlen_adjusted = int(np.median([b[0] for b in size_pos_bounded]))
            posA_adjusted = int(np.median([b[1] for b in size_pos_bounded]))
            posB_adjusted = int(np.median([b[2] for b in size_pos_bounded]))
            if svtype_m == 2:  # del
                posA = posA_adjusted
                posB = posB_adjusted
                svlen = svlen_adjusted
                er.preciseA = True
                er.preciseB = True
            else:
                # to_assemble = False
                er.preciseA = True
                er.preciseB = True
                posA = posA_adjusted
                posB = posB_adjusted
                svlen = svlen_adjusted

        else:
            svlen = int(np.median([sp[4] for sp in spanning_alignments]))
            posA = spanning_alignments[best_index][2]
            posB = spanning_alignments[best_index][3]
            er.preciseA = True
            er.preciseB = True

        ab = abs(posB - posA)
        if svlen > 0:
            jitter = ((posA_95 + posB_95) / 2) / svlen
        else:
            jitter = -1

        er.svtype = "DEL" if svtype_m == 2 else "INS"
        er.join_type = "3to5"
        er.chrA = chrom
        er.chrB = chrom
        er.posA = posA
        er.posB = posB
        er.cipos95A = posA_95
        er.cipos95B = posB_95
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

        count_attributes2(informative_reads, [], [i[5] for i in spanning_alignments], insert_ppf, generic_insertions, er, paired_end)

        er.contig = None
        er.contig_left_weight = 0
        er.contig_right_weight = 0
        er.contig2 = None
        er.contig2_left_weight = 0
        er.contig2_right_weight = 0
        as1 = None
        as2 = None
        ref_bases = 0
        if to_assemble:
            if er.preciseA:
                as1 = assembler.base_assemble(u_reads, er.posA, 500)
                if as1:
                    ref_bases += assign_contig_to_break(as1, er, "A", spanning_alignments)
            if er.preciseB:
                as2 = assembler.base_assemble(v_reads, er.posB, 500)
                if as2:
                    ref_bases += assign_contig_to_break(as2, er, "B", 0)
            if not as1 and len(generic_insertions) > 0:
                as1 = assembler.base_assemble([item.read_a for item in generic_insertions], er.posA, 500)
                if as1:
                    ref_bases += assign_contig_to_break(as1, er, "A", 0)
        er.linked = 0
        er.block_edge = 0
        er.ref_bases = ref_bases
        er.sqc = -1  # not calculated
        if er.svtype == "INS":
            cigar_index = spanning_alignments[best_index][6]
            start_ins = 0
            ct = best_align.cigartuples
            target_len = svlen
            # if er.preciseA:
            #     target_len = svlen
            # else:
            #     target_len = ct[cigar_index][1]
            for i in range(cigar_index):
                if ct[i][0] in {0, 1, 4, 7, 8}:
                    start_ins += ct[i][1]
            er.variant_seq = best_align.seq[start_ins:start_ins+target_len]
            er.ref_seq = best_align.seq[start_ins - 1]
        return er
    elif len(informative) > 0:
        return partition_single(informative, insert_size, insert_stdev, insert_ppf, spanning_alignments,
                                min_support, to_assemble, [], sites_info, paired_end)
    elif len(generic_insertions) > 0:
        if sites_info:
            site = sites_info[0]
        else:
            site = None
        # this is a bit confusing, generic insertions are used instead to make call, but bnd count is now 0
        return make_single_call([], insert_size, insert_stdev, insert_ppf, min_support, to_assemble,
                                spanning_alignments, 0, generic_insertions, site, paired_end)


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


cdef tuple break_ops(positions, precise, int limit, float median_pos):
    # Inspired by mosdepth algorithm +1 for start -1 for end using intervals where break site could occur
    # Use insert size to set a limit on where break site could occur
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


def value_closest_to_m(l):
    # cdef float svlen_m = np.mean(l)
    svlen_m = np.median(l)
    return min(l, key=lambda x:abs(x-svlen_m))


cdef void set_ins_seq(EventResult_t er):
    # if possible extract the insertion seq from alignment
    if not er.svtype == "INS":
        return
    cdef int svlen = er.svlen
    cdef int lc, lc2, start_i, end_i, start2_i, end2_i
    cdef bint done
    if er.svlen_precise:
        if er.join_type == "5to3" or er.join_type == "3to5":
            if er.contig and er.contig2:
                lc = len(er.contig)
                lc2 = len(er.contig2)
                if er.contig_ref_start < er.contig2_ref_end:
                    # use the aligned portion of contig as ins seq
                    # count back from right clip
                    done = False
                    cont2 = er.contig2
                    end2_i = lc2 - 1
                    while cont2[end2_i].islower():
                        end2_i -= 1
                        if end2_i < 0:
                            break
                    start2_i = end2_i + 1 - svlen
                    if start2_i >= 0:
                        seq = cont2[start2_i: end2_i]
                        er.variant_seq = seq
                        done = True
                    # try and use the left contig instead
                    if not done:
                        cont = er.contig
                        start_i = 0
                        while cont[start_i].islower():
                            start_i += 1
                            if start_i == lc:
                                break
                        end_i = start_i + svlen
                        if end_i < (lc - start_i):
                            seq = cont[start_i: end_i]
                            er.variant_seq = seq
                        else:
                            if start_i > 20:
                                er.left_ins_seq = cont[:start_i]
                            if lc2 - end2_i > 20:
                                er.right_ins_seq = cont2[-end2_i:]
                # else:
                #     happens rarely

            else:
                # try and use the soft clip portion of one, count from aligned base outwards
                if er.contig and (er.contig[0].isupper() or er.contig[-1].isupper()):  # skip both ends clipped
                    lc = len(er.contig)
                    if er.contig[0].isupper() and er.contig[-1].islower():  # fetch right
                        if er.contig_rc >= svlen:
                            er.variant_seq = er.contig[lc - er.contig_rc: lc - er.contig_rc + svlen]
                        else:
                            er.right_ins_seq = er.contig[lc - er.contig_rc:]
                    elif er.contig[-1].isupper() and er.contig[0].islower():   # fetch left
                        if er.contig_lc >= svlen:
                            er.variant_seq = er.contig[er.contig_lc - svlen: er.contig_lc]
                        else:
                            er.left_ins_seq = er.contig[:er.contig_lc]
                elif er.contig2 and (er.contig2[0].isupper() or er.contig2[-1].isupper()):
                    lc2 = len(er.contig2)
                    if er.contig2[0].isupper() and er.contig2[-1].islower():
                        if er.contig2_rc >= svlen:
                            er.variant_seq = er.contig2[lc2 - er.contig2_rc: lc2 - er.contig2_rc + svlen]
                        else:
                            er.right_ins_seq = er.contig2[lc2 - er.contig2_rc:]
                    elif er.contig2[-1].isupper() and er.contig2[0].islower():
                        if er.contig2_lc >= svlen:
                            er.variant_seq = er.contig2[er.contig2_lc - svlen: er.contig2_lc]
                        else:
                            er.left_ins_seq = er.contig2[:er.contig2_lc]


cdef infer_unmapped_insertion_break_point(int main_A_break, int cipos95A, int preciseA, int main_B_break, int cipos95B,
                                          int preciseB):
    # A or B break pay be in accurate, try and infer svlen and a more accurate break point(s)
    cdef int svlen = 0
    cdef int mid
    if not preciseA and not preciseB:
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
    cdef int limit = insert_size + insert_stdev
    cdef AlignmentItem i
    positionsA = [i.breakA for i in informative]
    positionsB = [i.breakB for i in informative]
    cdef float median_A = np.median(positionsA)
    cdef float median_B = np.median(positionsB)
    cdef int main_A_break = 0
    cdef int main_B_break = 0
    cdef int cipos95A, cipos95B
    cdef float svlen_m
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

    svlen_precise = 0
    if informative[0].chrA == informative[0].chrB:
        if svtype == "INS":  # use inferred query gap to call svlen
            lens = []
            inferred_lens = []
            for i in informative:
                if i.inferred_sv_len != -1:
                    if i.size_inferred == 1:
                        inferred_lens.append(i.inferred_sv_len)
                    else:
                        lens.append(i.inferred_sv_len)
            if len(lens) > 0:
                svlen = value_closest_to_m(lens)
                svlen_precise = 1
            elif len(inferred_lens) > 0:
                svlen = value_closest_to_m(inferred_lens)
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
                    svlen = value_closest_to_m(lens)
                    if main_svlen > 0 and (svlen / main_svlen) > 0.7:
                        svlen_precise = 1
                    else:
                        svlen = main_svlen
                else:
                    if len(inferred_lens) > 0:
                        svlen = value_closest_to_m(inferred_lens)
                    else:
                        svlen = main_svlen

                    if main_svlen > 0 and (svlen / main_svlen) <= 0.7:
                        svlen = main_svlen
        else:
            svlen = abs(main_B_break - main_A_break)
    else:
        svlen = -1
    if svlen == 0:
        svlen = -1
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
    cl = (4, 5)
    cdef int left_clipA = a_ct[0][1] if a_ct[0][0] in cl else 0
    cdef int right_clipA = a_ct[-1][1] if a_ct[-1][0] in cl else 0
    cdef int left_clipB = b_ct[0][1] if b_ct[0][0] in cl else 0
    cdef int right_clipB = b_ct[-1][1] if b_ct[-1][0] in cl else 0
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
    cdef int end = 0
    cdef int start = 0
    cdef int query_length = r.infer_read_length()  # Note, this also counts hard-clips
    end = query_length
    if r.cigartuples[0][0] == 4 or r.cigartuples[0][0] == 5:
        start += r.cigartuples[0][1]
    if r.cigartuples[-1][0] == 4 or r.cigartuples[-1][0] == 5:
        end -= r.cigartuples[-1][1]
    return start, end, query_length


cdef start_end_query_pair(r1, r2):
    cdef int query_length, s1, e1, s2, e2, start_temp, r1l, r2l
    # r1 and r2 might be on same read, if this is case infer the query position on the read
    s1, e1, r1l = query_start_end_from_cigartuples(r1)
    s2, e2, r2l = query_start_end_from_cigartuples(r2)
    if r1.flag & 64 == r2.flag & 64:  # same read
        if r2.flag & 16 != r1.flag & 16:  # different strand, count from end
            query_length = r1l  # Note, this also counts hard-clips
            start_temp = query_length - e2
            e2 = start_temp + e2 - s2
            s2 = start_temp
    return s1, e1, s2, e2, r1l, r2l


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
        er.contig_lc = as1["left_clips"]
        er.contig_rc = as1["right_clips"]
    if as2 is not None and "contig" in as2:
        er.contig2 = as2["contig"]
        rbases += as2["ref_bases"]
        er.contig2_ref_start = as2["ref_start"]
        er.contig2_ref_end = as2["ref_end"]
        er.contig2_left_weight = as2["left_weight"]
        er.contig2_right_weight = as2["right_weight"]
        er.contig2_lc = as2["left_clips"]
        er.contig2_rc = as2["right_clips"]
    er.ref_bases = rbases


cdef call_from_reads(u_reads, v_reads, int insert_size, int insert_stdev, float insert_ppf, int min_support,
                     int block_edge, int assemble, info, bint paired_end):
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
            a_qstart, a_qend, b_qstart, b_qend, a_len, b_len = start_end_query_pair(a, b)
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
                                   a.pos, a.reference_end,
                                   b.pos, b.reference_end,
                                   3 if a.flag & 16 == 0 else 5,
                                   3 if b.flag & 16 == 0 else 5,
                                   left_clip_a, right_clip_a,
                                   left_clip_b, right_clip_b,
                                   a_qstart, a_qend, b_qstart, b_qend, a_len, b_len,
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
            count_attributes2(u_reads, v_reads, [], insert_ppf, [], er, paired_end)
            if er.su < min_support:
                continue
            er.svlen_precise = 1  # prevent remapping
            assemble_partitioned_reads(er, u_reads, v_reads, block_edge, assemble)
            if info:
                info = [st for st in info if st.svtype == er.svtype]
                if len(info) > 0:
                    if len(info) == 1:
                        site = info[0]
                    else:
                        best_i = 0
                        best_sep = 1000000000
                        for si, site in enumerate(info):
                            sep = np.sqrt((site.start - er.posA)**2 + (site.end - er.posB)**2)
                            if sep < best_sep:
                                best_sep = sep
                                best_i = si
                        site = info[best_i]
                    er.site_info = site
            results.append(er)

    return results


cdef filter_single_partitions(u_reads, v_reads):
    # rare, but single reads with >2 alignments can have multiple alignments end up in one block. These should be processed as singles
    u_counts = defaultdict(list)
    v_counts = defaultdict(list)
    any_u_grouped = False
    any_v_grouped = False
    for cigar_info, a in u_reads:
        u_counts[(a.is_read1, a.qname)].append((cigar_info, a))
        if len(u_counts[(a.is_read1, a.qname)]) > 1:
            any_u_grouped = True
    for cigar_info, a in v_reads:
        v_counts[(a.is_read1, a.qname)].append((cigar_info, a))
        if len(v_counts[(a.is_read1, a.qname)]) > 1:
            any_v_grouped = True
    if not any_u_grouped and not any_v_grouped:
        return u_reads, v_reads, None, None
    single_u, single_v, actual_u, actual_v = [], [], [], []
    for k, v in u_counts.items():
        if len(v) == 1:
            actual_u += v
        else:
            single_u += v
    for k, v in v_counts.items():
        if len(v) == 1:
            actual_v += v
        else:
            single_v += v
    return actual_u, actual_v, single_u, single_v


cdef one_edge(u_reads_info, v_reads_info, int clip_length, int insert_size, int insert_stdev, float insert_ppf,
                   int min_support, int block_edge, int assemble, info, bint paired_end):
    spanning_alignments = []
    u_reads = []
    v_reads = []
    cdef int cigar_index, event_pos
    for cigar_info, a in u_reads_info:
        if not a.cigartuples:
            continue
        u_reads.append(a)
        cigar_index = cigar_info.cigar_index
        if 0 < cigar_index < len(a.cigartuples) - 1:  # Alignment spans SV
            event_pos = cigar_info.event_pos
            ci = a.cigartuples[cigar_index]
            spanning_alignments.append((ci[0],
                                        a.rname,
                                        event_pos,
                                        event_pos + 1 if ci[0] == 1 else event_pos + a.cigartuples[cigar_index][1] + 1,
                                        event_pos + a.cigartuples[cigar_index][1] + 1,
                                        a,
                                        cigar_index))

    for cigar_info, a in v_reads_info:
        if not a.cigartuples:
            continue
        v_reads.append(a)
        cigar_index = cigar_info.cigar_index
        if 0 < cigar_index < len(a.cigartuples) - 1:  # Alignment spans SV
            event_pos = cigar_info.event_pos
            ci = a.cigartuples[cigar_index]
            spanning_alignments.append((ci[0],
                                        a.rname,
                                        event_pos,
                                        event_pos + 1 if ci[0] == 1 else event_pos + a.cigartuples[cigar_index][1] + 1,
                                        event_pos + a.cigartuples[cigar_index][1],
                                        a,
                                        cigar_index))

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
        # choose representative alignment to use
        best_index = 0
        best_dist = 1e9
        for index in range(len(spanning_alignments)):
            dist = abs(spanning_alignments[index][2] - posA) + abs(spanning_alignments[index][3] - posB)
            if dist < best_dist:
                best_index = index
                best_dist = dist
                if dist == 0:
                    break

        best_align = spanning_alignments[best_index][5]
        svlen = spanning_alignments[best_index][4]
        posA = spanning_alignments[best_index][2]
        posB = spanning_alignments[best_index][3]
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
        er.jitter = jitter
        er.preciseA = True
        er.preciseB = True
        er.svlen = svlen
        er.query_overlap = 0
        u_reads = [i[5] for i in spanning_alignments]
        v_reads = []
        count_attributes2(u_reads, v_reads, [i[5] for i in spanning_alignments], insert_ppf, [], er, paired_end)
        if er.su < min_support:
            return []

        assemble_partitioned_reads(er, u_reads, v_reads, block_edge, assemble)
        if er.svtype == "INS":
            cigar_index = spanning_alignments[best_index][6]
            start_ins = 0
            ct = best_align.cigartuples
            for i in range(cigar_index):
                if ct[i][0] in {0, 1, 4, 7, 8}:
                    start_ins += ct[i][1]
            er.variant_seq = best_align.seq[start_ins:start_ins+svlen]
            er.ref_seq = best_align.seq[start_ins - 1]
        if info:
            info = [i for i in info if i.svtype == er.svtype]
            if len(info) > 0:
                site = info[0]
                if len(info) > 1:
                    best_i = 0
                    best_sep = 1000000000
                    for si, site in enumerate(info):
                        sep = np.sqrt((site.start - er.posA)**2 + (site.end - er.posB)**2)
                        if sep < best_sep:
                            best_sep = sep
                            best_i = si
                    site = info[best_i]
                er.site_info = site
        return [er]

    else:
        results = call_from_reads(u_reads, v_reads, insert_size, insert_stdev, insert_ppf, min_support, block_edge, assemble, info, paired_end)
        return results


def fpos_srt(x):
    return x[0].tell


cdef get_reads(infile, nodes_info, buffered_reads, n2n, bint add_to_buffer, sites_index):
    cdef int j, int_node, steps
    cdef uint64_t p
    cdef uint64_t v
    cdef AlignedSegment a
    aligns = []
    fpos = []
    site_info = []
    for int_node in nodes_info:
        if sites_index and int_node in sites_index:
            continue
        n = n2n[int_node]
        if int_node in buffered_reads:
            aligns.append((n, buffered_reads[int_node]))
            continue
        fpos.append((n, int_node))
    for node, int_node in sorted(fpos, key=fpos_srt):
        infile.seek(node.tell)
        try:
            a = next(infile)
        except StopIteration:
            return aligns
        v = xxhasher(bam_get_qname(a._delegate), len(a.qname), 42)
        if v == node.hash_name and a.flag == node.flag and a.pos == node.pos and a.rname == node.chrom:
            aligns.append((node, a))
            if add_to_buffer:
                buffered_reads[int_node] = a  # Add to buffer, then block nodes with multi-edges dont need collecting twice
            continue
        else:  # Try next few reads, find the read in the bgzf block?
            steps = 0
            while steps < 50:
                try:
                    a = next(infile)
                except StopIteration:
                    return aligns
                steps += 1
                v = xxhasher(bam_get_qname(a._delegate), len(a.qname), 42)
                if v == node.hash_name and a.flag == node.flag and a.pos == node.pos and a.rname == node.chrom:
                    aligns.append((node, a))
                    if add_to_buffer:
                        buffered_reads[int_node] = a
                    break

    return aligns


cdef list multi(data, bam, int insert_size, int insert_stdev, float insert_ppf, int clip_length, int min_support, int lower_bound_support,
                 int assemble_contigs, int max_single_size, info, bint paired_end, int length_extend, float divergence):

    # Sometimes partitions are not linked, happens when there is not much support between partitions
    # Then need to decide whether to call from a single partition
    n2n = data["n2n"]
    seen = set(range(len(data['parts'])))
    out_counts = defaultdict(int)  # The number of 'outward' links to other clusters
    cdef int buffered_reads = 0
    cdef bint add_to_buffer = 1
    cdef int int_node
    cdef unsigned long[:] partition
    events = []
    if info:
        sites_info = list(info.values())
    else:
        sites_info = []
    # u and v are the part ids, d[0] and d[1] are the lists of nodes for those parts
    for (u, v), d in data["s_between"].items():
        rd_u = get_reads(bam, d[0], data["reads"], data["n2n"], add_to_buffer, info)   # [(Nodeinfo, alignment)..]
        rd_v = get_reads(bam, d[1], data["reads"], data["n2n"], add_to_buffer, info)
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

        # finds reads that should be a single partition
        u_reads, v_reads, u_single, v_single = filter_single_partitions(rd_u, rd_v)
        if len(u_reads) > 0 and len(v_reads) > 0:
            events += one_edge(rd_u, rd_v, clip_length, insert_size, insert_stdev, insert_ppf, min_support, 1, assemble_contigs,
                               sites_info, paired_end)
        if u_single:
            res = single(u_single, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs,
                         sites_info, paired_end, length_extend, divergence)
            if res:
                if isinstance(res, EventResult):
                    events.append(res)
                else:
                    events += res
        if v_single:
            res = single(v_single, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs,
                         sites_info, paired_end, length_extend, divergence)
            if res:
                if isinstance(res, EventResult):
                    events.append(res)
                else:
                    events += res

    # Process any singles / unconnected blocks
    if seen:
        for part_key in seen:
            d = data["parts"][part_key]
            lb = lower_bound_support if len(sites_info) == 0 else 1
            if max_single_size > len(d) >= lb:
                # Call single block, only collect local reads to the block
                rds = get_reads(bam, d, data["reads"], data["n2n"], 0, info)
                if len(rds) < lower_bound_support or (len(sites_info) != 0 and len(rds) == 0):
                    continue
                res = single(rds, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs,
                             sites_info, paired_end, length_extend, divergence)
                if res:
                    if isinstance(res, EventResult):
                        events.append(res)
                    else:
                        events += res

    # Check for events within clustered nodes
    for k, d in data["s_within"].items():
        o_count = out_counts[k]
        i_counts = len(d)
        if i_counts > max_single_size:
            continue
        if o_count > 0 and i_counts > (2*min_support) and i_counts > o_count:
            rds = get_reads(bam, d, data["reads"], data["n2n"], 0, info)
            if len(rds) < lower_bound_support or (len(sites_info) != 0 and len(rds) == 0):
                    continue
            res = single(rds, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs,
                         sites_info, paired_end, length_extend, divergence)
            if res:
                if isinstance(res, EventResult):
                    events.append(res)
                else:
                    events += res
    return events


cpdef list call_from_block_model(bam, data, clip_length, insert_size, insert_stdev, insert_ppf, min_support, lower_bound_support,
                                 assemble_contigs, max_single_size, sites_index, bint paired_end, int length_extend, float divergence):
    n_parts = len(data["parts"])
    events = []
    if "info" in data:
        info = data["info"]
    else:
        info = None
    # next deal with info - need to filter these into the partitions, then deal with them in single / one_edge
    cdef EventResult_t e
    if n_parts >= 1:
        events += multi(data, bam, insert_size, insert_stdev, insert_ppf, clip_length, min_support, lower_bound_support,
                        assemble_contigs, max_single_size, info, paired_end, length_extend, divergence)
    elif n_parts == 0:
        if len(data["n2n"]) > max_single_size:
            return []
        rds = get_reads(bam, data["n2n"].keys(), data["reads"], data["n2n"], 0, info)
        if info:
            sites_info = list(info.values())
        else:
            sites_info = []
        if len(rds) < lower_bound_support or (len(sites_info) != 0 and len(rds) == 0):
            return []
        ev = single(rds, insert_size, insert_stdev, insert_ppf, clip_length, min_support, assemble_contigs, sites_info, paired_end, length_extend, divergence)
        if ev:
            if isinstance(ev, list):
                events += ev
            else:
                events.append(ev)
    events = [e for e in events if e and (e.svlen > 0 or e.svtype == "TRA")]
    for e in events:
        if e.svlen_precise:
            set_ins_seq(e)
    return events
