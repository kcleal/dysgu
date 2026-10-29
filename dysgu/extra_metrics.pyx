#cython: language_level=3, c_string_type=unicode, c_string_encoding=utf8
import logging
import multiprocessing
import time
import pysam
import numpy as np
cimport numpy as np
from dysgu.map_set_utils cimport EventResult, cigar_clip, clip_sizes, clip_sizes_hard
from dysgu.map_set_utils cimport TranscriptData
from dysgu.map_set_utils cimport map as ankerl_map
from dysgu.map_set_utils cimport set as ankerl_set
from dysgu.map_set_utils cimport hash as xxhasher
from dysgu.map_set_utils import merge_intervals
from dysgu.io_funcs import intersecter #, iitree
from superintervals import IntervalMap
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement
from libc.math cimport fabs as c_fabs, log as c_log, exp as c_exp
from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.vector cimport vector as cpp_vector
from dysgu.map_set_utils import echo
import array
from collections import defaultdict
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_cigar, bam_get_qname
import pandas as pd
pd.options.mode.chained_assignment = None
from sys import byteorder
import os
import resource
ctypedef EventResult EventResult_t


class BadClipCounter:

    def __init__(self, n_references, low_mem, temp_dir):
        self.low_mem = low_mem
        self.temp_dir = temp_dir
        self.n_refs = n_references
        if not self.low_mem:
            self.clip_pos_arr = [array.array("L", []) for i in range(n_references)]  # 32 bit unsigned int
        else:
            soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
            if soft + (n_references * 2) > soft:
                resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, soft + (n_references * 2)), hard))
            self.clip_pos_arr = [open(f"{self.temp_dir}/{i}.badclip.bin", "wb") for i in range(n_references)]

    def tidy(self):
        if self.low_mem:
            [os.remove(f"{self.temp_dir}/{i}.badclip.bin") for i in range(self.n_refs)]

    def add(self, chrom, position):
        if not self.low_mem:
            self.clip_pos_arr[chrom].append(position)
        else:
            self.clip_pos_arr[chrom].write(position.to_bytes(4, byteorder))

    def sort_arrays(self):
        if not self.low_mem:
            self.clip_pos_arr = [np.array(i, dtype=np.dtype("i")) for i in self.clip_pos_arr]
        else:
            [i.close() for i in self.clip_pos_arr]
            self.clip_pos_arr = [np.fromfile(f"{self.temp_dir}/{i}.badclip.bin", np.int32) for i in range(self.n_refs)]

        [i.sort() for i in self.clip_pos_arr]

    def count_near(self, int chrom, int start, int end):
        py_a = self.clip_pos_arr[chrom]
        cdef int len_a = len(py_a)
        if chrom < 0 or len_a == 0:
            return 0
        cdef int[:] a = py_a
        if start < 0:
            start = 0
        cdef int i = np.searchsorted(py_a, start)
        if i < 0:
            i = 0
        if i >= len_a:
            return 0
        # search forward and backwards
        cdef int count = 0
        cdef int p
        while True:
            if i >= len_a:
                break
            p = a[i]
            if p <= end:
                count += 1
                i += 1
            else:
                break
        return count


cdef float soft_clip_qual_corr(reads):
    """Function to compute the correlation of quality values between the soft-clipped portions of reads. values of
    1 imply no correlation whereas closer to 0 implies that quality values are correlated"""
    cdef const unsigned char[:] quals
    cdef int idx, x, i, j
    cdef ankerl_map[int, cpp_vector[int]] qq
    cdef int left_clip, right_clip
    cdef AlignedSegment r
    for r in reads:
        quals = r.query_qualities
        left_clip = 0
        right_clip = 0
        clip_sizes(r, &left_clip, &right_clip)

        if quals is None:
            continue
        if left_clip:
            idx = 0
            for x in range(r.pos - left_clip, r.pos):
                qq[x].push_back(quals[idx])
                idx += 1
        if right_clip:
            end_pos = r.reference_end + right_clip - 1
            for idx in range(-1, - right_clip, -1):
                qq[end_pos].push_back(quals[idx])
                end_pos -= 1

    if qq.empty():
        return -1

    cdef float all_z = 0
    cdef bint seen = False
    cdef cpp_vector[int] all_v
    cdef float sum_all_v = 0
    cdef cpp_vector[int]* second_ptr
    cdef int first
    cdef float mean_second, sum_second, sum_diff
    cdef ankerl_map[int, cpp_vector[int]].iterator qq_iter
    cdef int size;
    with nogil:
        qq_iter = qq.begin()
        while qq_iter != qq.end():
            second_ptr = &dereference(qq_iter).second
            size = second_ptr[0].size()
            if size > 1:
                sum_second = 0
                for i in range(size):
                    sum_second += second_ptr[0][i]
                mean_second = sum_second / float(size)
                for i in range(size):
                    all_z += c_fabs(second_ptr[0][i] - mean_second)
                sum_all_v += sum_second
                all_v.insert(all_v.end(), second_ptr[0].begin(), second_ptr[0].end())  # concat
                seen = True
            preincrement(qq_iter)
    if not seen:
        return -1

    cdef float mean_all_v = sum_all_v / float(all_v.size())
    cdef float z = 0
    for i in range(all_v.size()):
        z += c_fabs(all_v[i] - mean_all_v)
    if z == 0:
        return -1

    return all_z / z


cdef int CMATCH = 0
cdef int CINS = 1
cdef int CDEL = 2
cdef int CREF_SKIP = 3
cdef int CSOFT_CLIP = 4
cdef int CHARD_CLIP = 5
cdef int CEQUAL = 7
cdef int CDIFF = 8

cdef int BAM_CIGAR_SHIFT = 4
cdef int BAM_CIGAR_MASK  = 15 #0xf


cdef inline int i_min(int x, int y) noexcept nogil:
    return x if x < y else y


cdef inline int i_max(int x, int y) noexcept nogil:
    return x if x > y else y


cdef inline int i_abs(int x) noexcept nogil:
    return x if x >= 0 else -x


cdef extern from "si_intmap.h":
    cdef cppclass SiIntMap:
        SiIntMap() except +
        void add(int start, int end, int value) nogil
        void build() except + nogil
        void clear() nogil
        void search_values(int start, int end, cpp_vector[int]& found) nogil const


cdef inline uint64_t aln_qhash(AlignedSegment aln) noexcept:
    return xxhasher(bam_get_qname(aln._delegate), aln._delegate.core.l_qname, 42)


cdef inline uint64_t support_key(int event_idx, uint64_t qhash) noexcept nogil:
    return qhash ^ ((<uint64_t>(event_idx + 1)) * (<uint64_t>11400714819323198485))


cdef inline bint add_support_hash(ankerl_set[uint64_t]& values, int event_idx, uint64_t qhash) noexcept nogil:
    cdef uint64_t key = support_key(event_idx, qhash)
    if values.find(key) == values.end():
        values.insert(key)
        return True
    return False


def _bam_chrom_label(chrom, references, chrom_to_tid):
    if chrom in chrom_to_tid:
        return chrom
    if isinstance(chrom, bytes):
        chrom = chrom.decode()
        if chrom in chrom_to_tid:
            return chrom
    elif not isinstance(chrom, str):
        try:
            tid = int(chrom)
        except (TypeError, ValueError):
            tid = -1
        if 0 <= tid < len(references):
            return references[tid]
    return chrom


cdef inline bint cigar_op_consumes_ref(int op) noexcept nogil:
    return op == CMATCH or op == CDEL or op == CREF_SKIP or op == CEQUAL or op == CDIFF


def recount_support_from_bam(events, bam, int mapq_threshold=0, int read_length=150, int insert_median=-1,
                             int insert_stdev=-1, bint paired_end=True, int max_fetch_window=2000,
                             int breakpoint_window=50, gaps_file=""):
    if not events:
        return events

    t0 = time.time()
    try:
        bam.check_index()
    except (AttributeError, ValueError):
        logging.info("Skipping SKIP support recount: --ibam is not indexed")
        return events

    cdef int n_events = len(events)
    cdef int n_skip = 0
    cdef int updated = 0
    cdef int alignments_seen = 0
    cdef int reads_seen = 0
    cdef int pad = 50
    cdef int idx, ta, tb, start, end, chrom_len, map_start, map_end
    cdef int tid, roi_start, roi_end, roi_s, roi_e, ref_pos, i, op, length, gap_start, gap_end, event_idx
    cdef int event_start, event_end, event_len, gap_len, overlap_start, overlap_end, overlap
    cdef int ann_filtered = 0
    cdef int best_idx, best_support, ev_support
    cdef double ann_score = 0.0, sv_score = 0.0, best_score
    cdef bint have_gaps = False
    cdef str chrom_b
    cdef TranscriptData transcript_gaps
    cdef uint32_t cigar_l, cigar_value
    cdef uint32_t *cigar_p
    cdef uint64_t qhash
    cdef EventResult_t event
    cdef EventResult_t ev
    cdef AlignedSegment aln
    cdef ankerl_set[uint64_t] seen
    cdef object interval_map
    cdef object candidates
    wr_counts = [0] * n_events
    references = bam.references
    ref_lengths = bam.lengths
    chrom_to_tid = {chrom: tid for tid, chrom in enumerate(references)}
    skip_maps = [None] * len(references)
    skip_tids = [-1] * n_events
    skip_starts = [0] * n_events
    skip_ends = [0] * n_events
    map_records_by_tid = defaultdict(list)
    roi_by_tid = defaultdict(list)
    rois = []

    if gaps_file:
        transcript_gaps = TranscriptData()
        transcript_gaps.readBed(gaps_file)
        have_gaps = True

    for idx, event in enumerate(events):
        if event.svtype != "SKIP":
            continue
        chrom_a = _bam_chrom_label(event.chrA, references, chrom_to_tid)
        chrom_b = _bam_chrom_label(event.chrB, references, chrom_to_tid)
        if chrom_a != event.chrA:
            event.chrA = chrom_a
        if chrom_b != event.chrB:
            event.chrB = chrom_b
        ta = chrom_to_tid.get(chrom_a, -1)
        tb = chrom_to_tid.get(chrom_b, -1)
        if ta < 0 or ta != tb:
            continue
        start = i_min(event.posA, event.posB)
        end = i_max(event.posA, event.posB)
        if end <= start:
            continue
        chrom_len = ref_lengths[ta]
        map_start = max(0, event.posA - pad)
        map_end = min(chrom_len, event.posA + pad + 1)
        if map_end <= map_start:
            continue
        map_records_by_tid[ta].append((map_start, map_end, idx))
        roi_by_tid[ta].append((map_start, map_end))
        skip_tids[idx] = ta
        skip_starts[idx] = start
        skip_ends[idx] = end
        n_skip += 1

    for ta, records in map_records_by_tid.items():
        interval_map = IntervalMap()
        for map_start, map_end, idx in records:
            interval_map.add(map_start, map_end, idx)
        interval_map.build()
        skip_maps[ta] = interval_map

    for tid, intervals in roi_by_tid.items():
        intervals.sort()
        roi_start, roi_end = intervals[0]
        for roi_s, roi_e in intervals[1:]:
            if roi_s <= roi_end:
                if roi_e > roi_end:
                    roi_end = roi_e
            else:
                rois.append((tid, roi_start, roi_end))
                roi_start, roi_end = roi_s, roi_e
        rois.append((tid, roi_start, roi_end))
    rois.sort()

    if not rois:
        logging.info("Skipping SKIP support recount: no SKIP breakpoints found in --ibam references")
        return events

    for tid, roi_start, roi_end in rois:
        interval_map = skip_maps[tid]
        if interval_map is None:
            continue
        if have_gaps:
            chrom_b = references[tid]
        for aln in bam.fetch(references[tid], roi_start, roi_end):
            if aln.is_unmapped or aln.is_duplicate:
                continue
            alignments_seen += 1
            if aln._delegate.core.tid != tid:
                continue
            cigar_l = aln._delegate.core.n_cigar
            if cigar_l == 0:
                continue
            cigar_p = bam_get_cigar(aln._delegate)
            qhash = aln_qhash(aln)
            reads_seen += 1
            ref_pos = aln._delegate.core.pos
            for i in range(cigar_l):
                cigar_value = cigar_p[i]
                op = <int>cigar_value & BAM_CIGAR_MASK
                length = <int>cigar_value >> BAM_CIGAR_SHIFT
                if op == CREF_SKIP:
                    gap_start = ref_pos
                    gap_end = ref_pos + length
                    candidates = interval_map.search_values(gap_start, gap_start + 1)
                    if candidates:
                        if have_gaps:
                            ann_score = transcript_gaps.bestRefSkipOverlap(chrom_b, gap_start, gap_end)
                        # Credit only the candidate this skip matches best; ties go to the
                        # event with more independent support so zero-support duplicates
                        # are not resurrected by the recount
                        best_idx = -1
                        best_score = -1.0
                        best_support = -1
                        for event_idx in candidates:
                            if skip_tids[event_idx] != tid:
                                continue
                            event_start = skip_starts[event_idx]
                            event_end = skip_ends[event_idx]
                            event_len = event_end - event_start
                            gap_len = gap_end - gap_start
                            if event_len <= 0 or gap_len <= 0:
                                continue
                            overlap_start = i_max(gap_start, event_start)
                            overlap_end = i_min(gap_end, event_end)
                            overlap = overlap_end - overlap_start
                            if overlap <= 0:
                                continue
                            if overlap * 100 < event_len * 96 or overlap * 100 < gap_len * 96:
                                continue
                            sv_score = <double>overlap / <double>(gap_len if gap_len > event_len else event_len)
                            if have_gaps:
                                # Count support only when the SV explains this skip better than any known intron
                                if sv_score <= ann_score:
                                    ann_filtered += 1
                                    continue
                            ev = events[event_idx]
                            ev_support = ev.pe + ev.supp + ev.sc
                            if sv_score > best_score or (sv_score == best_score and ev_support > best_support):
                                best_score = sv_score
                                best_idx = event_idx
                                best_support = ev_support
                        if best_idx != -1:
                            if add_support_hash(seen, best_idx, qhash):
                                wr_counts[best_idx] += 1
                    ref_pos = gap_end
                elif cigar_op_consumes_ref(op):
                    ref_pos += length

    if alignments_seen == 0:
        # E.g. the --ibam index is missing or stale: region fetches silently return nothing
        # while sequential iteration still works. Do not clobber graph-derived support.
        logging.warning("SKIP support recount visited 0 alignments across {} ROIs; --ibam index may be missing or stale, keeping graph-derived support".format(len(rois)))
        return events

    for event_idx, event in enumerate(events):
        if event.svtype != "SKIP" or skip_tids[event_idx] < 0:
            continue
        if event.spanning != wr_counts[event_idx]:
            updated += 1
        event.spanning = wr_counts[event_idx]
        event.su = event.pe + event.supp + (2 * event.spanning) + event.bnd

    logging.info("SKIP recount: {} SKIP events ({} annotation-explained)".format(n_skip, ann_filtered))
    logging.debug("SKIP recount detail: {:.2f}s, {} WR changed, {} ROIs, {} reads ({} alignments visited)".format(
        time.time() - t0, updated, len(rois), reads_seen, alignments_seen))
    return events


cdef int GT_MAX_SPAN_ROI = 50000
cdef int GT_JOIN_SLOP = 250   # SA split junction tolerance
cdef int GT_FACE_SLOP = 200   # pair facing-end to breakpoint tolerance
cdef int GT_CLIP_MIN = 8      # genotyping clip floor; repeat insertions often clip shorter
                              # than the call-time clip_l, and the (0.75, 0.10) weight
                              # tolerates the higher false-clip rate

# Soft-likelihood genotype model. Each fragment carries an asymmetric error
# model (a, r): a = P(observation | alt allele), r = P(observation | ref
# allele), independent per evidence class (a clipped read has r = 0.10 even
# though a = 0.75). For genotype g with alt fraction pi_g a fragment adds
# log(a*pi_g + r*(1-pi_g)) to LL(g); pi = 0.005 (ref), 0.5 (het), 0.99
# (hom-alt). Under het every fragment contributes exactly log((a+r)/2), which
# is ~log(0.5) for decisive evidence, so LL_het is approximated from the
# fragment count and only LL_ref / LL_hom need to be accumulated per event.
PI_REF = 0.005
PI_HET = 0.5
PI_HOM = 0.99


cdef enum:
    SVT_DEL = 0
    SVT_INS = 1
    SVT_DUP = 2
    SVT_SKIP = 3
    SVT_OTHER = 4


cdef inline int _c_svt_code(str svtype) noexcept:
    if svtype == "DEL":
        return SVT_DEL
    if svtype == "INS":
        return SVT_INS
    if svtype == "DUP":
        return SVT_DUP
    if svtype == "SKIP":
        return SVT_SKIP
    return SVT_OTHER


cdef inline bint _c_len_match(int got, int want, int abs_tol, double frac) noexcept nogil:
    if want <= 0:
        return False
    return i_abs(got - want) <= i_max(abs_tol, <int>(frac * want))


cdef inline bint _gt_within(int v, int lo, int hi) noexcept nogil:
    return v >= lo and v <= hi


cdef inline double _gt_tlen_lr(int tlen, int mu, int sigma, int svlen, int direction) noexcept nogil:
    # Gaussian likelihood ratio for a spanning proper pair:
    # p = P(tlen | alt) / (P(tlen | alt) + P(tlen | ref)) with alt mean
    # mu + direction*svlen (DEL/INS lengthen the fragment, DUP shortens it).
    # Returned p is capped to [0.05, 0.98]; -1 if not applicable.
    # The log-odds are scaled by min(1, svlen / 2*sigma): below ~2 insert
    # stdev the alt and ref tlen distributions overlap so heavily that the
    # raw LR is per-fragment noise (the 0.98 cap would turn a sub-sigma
    # fluctuation into near-certain alt evidence, and het loci - where alt
    # and ref fragments straddle the midpoint - systematically read as hom).
    # Shrinking toward p = 0.5 keeps the vote direction but bounds its
    # weight, so uninformative fragments contribute ~log(0.5) under het.
    cdef double d_ref, d_alt, diff, p, w
    if tlen <= 0 or svlen <= 0 or mu <= 0 or sigma <= 0:
        return -1.0
    d_ref = (tlen - mu) / <double> sigma
    d_alt = (tlen - (mu + direction * svlen)) / <double> sigma
    # log(p / (1-p)) = 0.5 * (d_ref^2 - d_alt^2); evaluated in a form that
    # cannot overflow for tail observations
    diff = 0.5 * (d_ref * d_ref - d_alt * d_alt)
    w = svlen / (2.0 * <double> sigma)
    if w > 1.0:
        w = 1.0
    diff *= w
    if diff >= 0:
        p = 1.0 / (1.0 + c_exp(-diff))
    else:
        p = c_exp(diff) / (1.0 + c_exp(diff))
    if p < 0.05:
        p = 0.05
    elif p > 0.98:
        p = 0.98
    return p


cdef inline bint _c_split_match(const int* c1s, const int* p1s, const int* c2s, const int* p2s, const int* qgs,
                                int n_splits, int ta, int tb, int posA, int posB, int svlen, int svt, int ltol) noexcept nogil:
    # structural test for split-join alt evidence, mirroring how the event was
    # called: same-chrom joins must separate the two breakpoints by ~svlen in
    # reference (DEL/DUP/INV), insertion joins must separate the two alignments
    # by ~svlen in query coordinates. Without this, spurious short-range
    # supplementary chunks near a point event vote as strong alt evidence.
    cdef int i, p1, p2, qgap, c1, c2
    for i in range(n_splits):
        c1 = c1s[i]; p1 = p1s[i]; c2 = c2s[i]; p2 = p2s[i]; qgap = qgs[i]
        if c1 == ta and c2 == tb and i_abs(p1 - posA) <= GT_JOIN_SLOP and i_abs(p2 - posB) <= GT_JOIN_SLOP:
            if ta != tb:
                return True
            if svt == SVT_INS:
                if _c_len_match(qgap, svlen, ltol, 0.25):
                    return True
            elif _c_len_match(i_abs(p2 - p1), svlen, ltol, 0.15):
                return True
        if c1 == tb and c2 == ta and i_abs(p1 - posB) <= GT_JOIN_SLOP and i_abs(p2 - posA) <= GT_JOIN_SLOP:
            if ta != tb:
                return True
            if svt == SVT_INS:
                if _c_len_match(qgap, svlen, ltol, 0.25):
                    return True
            elif _c_len_match(i_abs(p2 - p1), svlen, ltol, 0.15):
                return True
    return False


cdef inline int _c_facing_end(int mtid, int tid, int mpos, int start, int aend, int bpos) noexcept nogil:
    # breakpoint-facing end of the read, observed coordinates only
    if mtid == tid:
        return aend if mpos > start else start
    return aend if i_abs(aend - bpos) <= i_abs(start - bpos) else start


cdef inline void _c_alt_weight(int ta, int tb, int posA, int posB, int start, int end, int span, int svlen, int svt,
                               int f_tid, int f_start, int f_aend, int f_ext_lo, int f_ext_hi,
                               int f_tlen, bint f_proper, bint f_paired, int f_clip_lo, int f_clip_hi,
                               int f_mtid, int f_mpos,
                               const int* op_v, const int* ln_v, const int* opos_v, int n_ops,
                               const int* c1s, const int* p1s, const int* c2s, const int* p2s, const int* qgs, int n_splits,
                               bint paired_end, int max_dist, int insert_median, int insert_stdev, int min_sv_size, int clip_l, int cslop,
                               double* out_a, double* out_r) noexcept nogil:
    """Alt-allele evidence for a fragment against an event as an (a, r) pair:
    a = P(observation | alt allele), r = P(observation | ref allele); (-1, -1)
    for no alt support. The asymmetric error model matters: a clipped read has
    r = 0.10 (clips happen by error) even though a = 0.75. Signature classes
    mirror the graph.pyx construct_graph decision tree so genotype-time evidence
    matches call-time evidence: CIGAR I/D/N, SA split joins (AlignmentsSA),
    breakend clips on proper pairs, and discordant pairs.
    CIGAR codes: M0 I1 D2 N3 S4 H5."""
    cdef int ltol = 10 if svlen < 150 else 50
    cdef double best = -1.0
    cdef double sim, p, v
    cdef int i, op, ln, opos, oend, overlap, direction
    cdef bint normal
    for i in range(n_ops):
        op = op_v[i]; ln = ln_v[i]; opos = opos_v[i]
        oend = opos + ln if (op == CDEL or op == CREF_SKIP) else opos
        if svt == SVT_DEL and op == CDEL:
            if _gt_within(opos, start - cslop, end + cslop) and _c_len_match(ln, svlen, ltol, 0.15):
                sim = 1.0 - c_fabs(<double>(ln - svlen)) / (svlen + 1.0)
                if sim < 0.0:
                    sim = 0.0
                v = 0.55 + 0.35 * sim
                if v > best:
                    best = v
        elif svt == SVT_INS and op == CINS:
            if not _gt_within(opos, start - cslop, end + cslop):
                continue
            if _c_len_match(ln, svlen, ltol, 0.25):
                sim = 1.0 - c_fabs(<double>(ln - svlen)) / (svlen + 1.0)
                if sim < 0.0:
                    sim = 0.0
                v = 0.55 + 0.35 * sim
                if v > best:
                    best = v
            elif 20 <= ln and ln < svlen - ltol:
                # partial insertion: the read captures only the start of the
                # inserted sequence (the rest runs off the read end)
                v = 0.5 + 0.3 * (<double> ln / (svlen + 1.0))
                if v > best:
                    best = v
        elif svt == SVT_SKIP and op == CREF_SKIP:
            if ln <= 0:
                continue
            overlap = i_min(oend, end) - i_max(opos, start)
            if overlap <= 0 or overlap * 100 < span * 96 or overlap * 100 < ln * 96:
                continue
            sim = <double>overlap / (ln if ln > span else span)
            v = 0.55 + 0.35 * sim
            if v > best:
                best = v
    if best >= 0:
        # CIGAR signature at the junction: specific evidence
        out_a[0] = best; out_r[0] = 0.05
        return
    # SA split joins (dysgu AlignmentsSA), either orientation of (posA, posB).
    # Long-read alignments yield more spurious supplementary joins, so the
    # false-positive rate r is higher when not paired-end.
    if n_splits and _c_split_match(c1s, p1s, c2s, p2s, qgs, n_splits, ta, tb, posA, posB, svlen, svt, ltol):
        out_a[0] = 0.95; out_r[0] = 0.02 if paired_end else 0.08
        return
    if not f_paired:
        out_a[0] = -1.0; out_r[0] = -1.0
        return
    if clip_l > 0 and f_proper and f_tlen < max_dist and f_mtid == f_tid:
        # proper pairs are evidence only when clipped at a breakpoint (graph.pyx BREAKEND path)
        if f_clip_lo >= GT_CLIP_MIN and i_abs(f_start - posA) <= cslop:
            out_a[0] = 0.75; out_r[0] = 0.10
            return
        if f_clip_hi >= GT_CLIP_MIN and i_abs(f_aend - posB) <= cslop:
            out_a[0] = 0.75; out_r[0] = 0.10
            return
        if f_clip_lo >= GT_CLIP_MIN and i_abs(f_start - posB) <= cslop and ta == tb:
            out_a[0] = 0.75; out_r[0] = 0.10
            return
        if f_clip_hi >= GT_CLIP_MIN and i_abs(f_aend - posA) <= cslop and ta == tb:
            out_a[0] = 0.75; out_r[0] = 0.10
            return
    # allele assignment for spanning proper pairs: Gaussian template-length
    # likelihood ratio (a genotype question; graph.pyx only asks whether the
    # pair is anomalous enough to be evidence at all). Uninformative pairs get
    # p ~ 0.5 and contribute near-identical constants to every state. Only
    # DEL/DUP shift reference-coordinate template length; an INS in the
    # unsequenced gap is invisible to tlen, so INS pairs get no tlen vote.
    if (svt == SVT_DEL or svt == SVT_DUP) and f_proper and f_mtid == f_tid and \
            f_ext_lo <= start - 10 and f_ext_hi >= end + 10:
        direction = 1 if svt == SVT_DEL else -1
        p = _gt_tlen_lr(f_tlen, insert_median, insert_stdev, svlen, direction)
        if p < 0.5:
            out_a[0] = -1.0; out_r[0] = -1.0  # ref-leaning pair: let the ref path claim it
            return
        out_a[0] = p; out_r[0] = 1.0 - p
        return
    normal = f_proper and f_tlen < max_dist and f_mtid == f_tid
    if normal:
        out_a[0] = -1.0; out_r[0] = -1.0
        return
    # discordant pair: facing end near one breakpoint, mate near the other.
    # Only meaningful when the event shifts template length by at least ~2
    # insert stdevs (or is trans-chrom); for small same-chrom events the
    # mate-proximity test is satisfied vacuously by any nearby pair.
    if ta != tb or svlen >= 2 * insert_stdev:
        if f_tid == ta and f_mtid == tb and i_abs(f_mpos - posB) <= max_dist:
            if i_abs(_c_facing_end(f_mtid, f_tid, f_mpos, f_start, f_aend, posA) - posA) <= GT_FACE_SLOP:
                out_a[0] = 0.85; out_r[0] = 0.03
                return
        if f_tid == tb and f_mtid == ta and i_abs(f_mpos - posA) <= max_dist:
            if i_abs(_c_facing_end(f_mtid, f_tid, f_mpos, f_start, f_aend, posB) - posB) <= GT_FACE_SLOP:
                out_a[0] = 0.85; out_r[0] = 0.03
                return
    out_a[0] = -1.0; out_r[0] = -1.0


cdef inline void _c_ref_weight(int ta, int tb, int posA, int posB, int start, int end, int span, int svlen, int svt,
                               int f_start, int f_aend, int f_ext_lo, int f_ext_hi, int f_tlen,
                               bint f_proper, int f_mtid, int f_tid,
                               const int* op_v, const int* ln_v, const int* opos_v, int n_ops,
                               const int* c1s, const int* p1s, const int* c2s, const int* p2s, const int* qgs, int n_splits,
                               bint paired_end, int max_dist, int insert_median, int insert_stdev, int min_sv_size, int clip_l,
                               double* out_a, double* out_r) noexcept nogil:
    """Soft ref-allele weight for a fragment against an event, or -1 for no ref
    support. Two evidence classes, both requiring no extractable signature:
    (1) breakpoint-window cover - the fragment (pair) or the read itself spans
    posA/posB +/- 100 bp cleanly; this is the only feasible ref evidence when
    the event is larger than the fragment reach (e.g. het DELs >= insert size);
    (2) whole-event cover - pair span or read covers [start, end], weighted by
    the Gaussian template-length likelihood ratio for DEL/INS/DUP."""
    cdef bint same_chrom = ta == tb
    cdef int win, reach, rwin
    cdef bint cover_a, cover_b, window_cover, read_win_a, read_win_b, read_window_cover
    cdef bint span_cover, pair_cover, read_span_cover, read_cover
    cdef int check_lo = 0, check_hi = 0
    cdef int op, ln, opos, oend, direction
    cdef int i
    cdef double p_alt
    # breakpoint-window cover by the full fragment (pair) or by the read alone.
    # The window must fit the fragment's actual reach: a 200 bp window can
    # never be covered by a 148 bp read, which previously left INS events with
    # no read-level ref evidence at all.
    reach = f_ext_hi - f_ext_lo
    win = 100 if reach >= 210 else reach // 2 - 5
    if win >= 15:
        cover_a = f_ext_lo <= posA - win and f_ext_hi >= posA + win
        cover_b = same_chrom and f_ext_lo <= posB - win and f_ext_hi >= posB + win
    else:
        cover_a = cover_b = False
    window_cover = cover_a or cover_b
    reach = f_aend - f_start
    rwin = 100 if reach >= 210 else reach // 2 - 5
    if rwin >= 15:
        read_win_a = f_start <= posA - rwin and f_aend >= posA + rwin
        read_win_b = same_chrom and f_start <= posB - rwin and f_aend >= posB + rwin
    else:
        read_win_a = read_win_b = False
    read_window_cover = read_win_a or read_win_b
    span_cover = same_chrom and f_ext_lo <= start - 10 and f_ext_hi >= end + 10
    if span_cover:
        pair_cover = True
        read_span_cover = f_start <= start - 10 and f_aend >= end + 10
        check_lo = start; check_hi = end
    elif window_cover or read_window_cover:
        pair_cover = window_cover
        read_span_cover = False
        if cover_a:
            win = 100 if (f_ext_hi - f_ext_lo) >= 210 else (f_ext_hi - f_ext_lo) // 2 - 5
            check_lo = posA - win; check_hi = posA + win
        elif read_win_a:
            check_lo = posA - rwin; check_hi = posA + rwin
        elif cover_b:
            win = 100 if (f_ext_hi - f_ext_lo) >= 210 else (f_ext_hi - f_ext_lo) // 2 - 5
            check_lo = posB - win; check_hi = posB + win
        else:
            check_lo = posB - rwin; check_hi = posB + rwin
    else:
        out_a[0] = -1.0; out_r[0] = -1.0
        return
    read_cover = read_span_cover or read_window_cover
    if n_splits:
        if paired_end:
            out_a[0] = -1.0; out_r[0] = -1.0
            return
        # noisy long reads carry spurious supplementary alignments; only splits
        # that would count as alt evidence for this event disqualify an
        # otherwise clean ref claim
        if _c_split_match(c1s, p1s, c2s, p2s, qgs, n_splits, ta, tb, posA, posB, svlen, svt,
                          10 if svlen < 150 else 50):
            out_a[0] = -1.0; out_r[0] = -1.0
            return
    # a read claiming to cover the breakpoint window must not carry a matching
    # indel there; partial insertions (read captures only the start of the
    # inserted sequence) are alt evidence, so for INS events any substantial
    # insertion in the window disqualifies the ref claim. Long-read error
    # profiles produce noise indels well above min_size, so the deletion floor
    # is raised when not paired-end.
    ins_contam = 20 if svt == SVT_INS else min_sv_size
    del_contam = min_sv_size if paired_end else max(min_sv_size, 100)
    for i in range(n_ops):
        op = op_v[i]; ln = ln_v[i]; opos = opos_v[i]
        oend = opos + ln if (op == CDEL or op == CREF_SKIP) else opos
        if (op == CDEL or op == CREF_SKIP) and ln >= del_contam and oend > check_lo and opos < check_hi:
            out_a[0] = -1.0; out_r[0] = -1.0
            return
        if op == CINS and ln >= ins_contam and _gt_within(opos, check_lo, check_hi):
            out_a[0] = -1.0; out_r[0] = -1.0
            return
        if clip_l > 0 and (op == CSOFT_CLIP or op == CHARD_CLIP) and ln >= clip_l and \
                (i_abs(opos - posA) <= 25 or (same_chrom and i_abs(opos - posB) <= 25)):
            out_a[0] = -1.0; out_r[0] = -1.0
            return
    if not paired_end:
        if read_cover:
            out_a[0] = 0.02; out_r[0] = 0.95
        else:
            out_a[0] = -1.0; out_r[0] = -1.0
        return
    # PE: read-level evidence needs no insert statistics; pair-level is gated
    if read_cover:
        out_a[0] = 0.02; out_r[0] = 0.95
        return
    if not pair_cover:
        out_a[0] = -1.0; out_r[0] = -1.0
        return
    if not (f_proper and f_tlen < max_dist and f_mtid == f_tid):
        out_a[0] = -1.0; out_r[0] = -1.0
        return
    if (svt == SVT_DEL or svt == SVT_DUP) and svlen > 0:
        direction = 1 if svt == SVT_DEL else -1
        p_alt = _gt_tlen_lr(f_tlen, insert_median, insert_stdev, svlen, direction)
        if p_alt > 0.5:
            out_a[0] = -1.0; out_r[0] = -1.0  # alt-leaning pair is not ref evidence
            return
        out_a[0] = p_alt; out_r[0] = 1.0 - p_alt
        return
    if svt == SVT_INS:
        # an insertion does not shift reference-coordinate tlen, so flanking
        # proper pairs are uninformative; read-level window cover is the only
        # pair-adjacent ref evidence and it was handled above
        out_a[0] = -1.0; out_r[0] = -1.0
        return
    out_a[0] = 0.05; out_r[0] = 0.9

def genotype_counts_from_bam(events, bam, int mapq_threshold=0, int min_sv_size=100, int max_dist=1000,
                             int insert_median=300, int insert_stdev=150, int clip_l=15, bint paired_end=True,
                             int procs=1, str bam_path="", ref_genome=None):
    """Re-count alt/ref supporting fragments per event for genotyping, using one batched
    scan over merged regions of interest (same architecture as recount_support_from_bam).
    Fragment classification mirrors the graph.pyx construct_graph definitions so
    genotype-time evidence is consistent with call-time evidence. No read_length
    anywhere: PE pair spans come from |tlen|, long reads are fetched by overlap.
    Sets event.gt_alt / event.gt_ref (fragment counts) and event.gt_ll_ref /
    event.gt_ll_hom (soft log-likelihood sums; LL_het is analytic, see module
    constants). Events with no alt evidence leave everything unset so the caller
    falls back to the default genotyping heuristics. Internal only, not serialized
    to VCF."""
    if not events:
        return events
    from dysgu.graph import AlignmentsSA  # lazy: graph.pyx imports extra_metrics at module level
    t0 = time.time()
    try:
        bam.check_index()
    except (AttributeError, ValueError):
        logging.info("Skipping fragment genotyping: --ibam is not indexed")
        return events

    n_events = len(events)
    references = bam.references
    ref_lengths = bam.lengths
    chrom_to_tid = {chrom: tid for tid, chrom in enumerate(references)}

    pad = max_dist if paired_end else 250
    cdef int cslop = 50

    cdef int start, end, aend, pos, op, ln, mtid, mpos, tlen, ext_lo, ext_hi, lo, i, vote, tid_c, k, run_max, c_alt, c_ref
    cdef int cigar_value
    cdef long long isize
    cdef bint proper, paired
    cdef uint32_t cigar_l
    cdef uint32_t *cigar_p
    cdef uint64_t qhash
    cdef AlignedSegment aln
    cdef cpp_vector[int] v_op, v_oln, v_opos
    cdef cpp_vector[int] v_cand, v_c1, v_p1, v_c2, v_p2, v_qg
    ev_info = []
    rois_by_tid = defaultdict(list)
    map_records_by_tid = defaultdict(list)
    cdef np.ndarray[np.int32_t, ndim=1] ev_ta = np.full(n_events, -1, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_tb = np.full(n_events, -1, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_posA = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_posB = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_start = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_end = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_span = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_svlen = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ev_support = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.uint8_t, ndim=1] ev_svt = np.zeros(n_events, dtype=np.uint8)
    for idx, event in enumerate(events):
        chrom_a = _bam_chrom_label(event.chrA, references, chrom_to_tid)
        ta = chrom_to_tid.get(chrom_a, -1)
        if ta < 0:
            ev_info.append(None)
            continue
        chrom_b = _bam_chrom_label(event.chrB, references, chrom_to_tid)
        tb = chrom_to_tid.get(chrom_b, -1)
        posA = event.posA
        posB = event.posB
        if posA <= posB:
            start = posA
            end = posB
        else:
            start = posB
            end = posA
        svlen = abs(event.svlen) if event.svlen and event.svlen > 0 else end - start
        ev_info.append({"svtype": event.svtype, "ta": ta, "tb": tb, "posA": posA, "posB": posB,
                        "start": start, "end": end, "span": end - start, "svlen": svlen,
                        "support": event.pe + event.supp + event.sc})
        ev_ta[idx] = ta
        ev_tb[idx] = tb
        ev_posA[idx] = posA
        ev_posB[idx] = posB
        ev_start[idx] = start
        ev_end[idx] = end
        ev_span[idx] = end - start
        ev_svlen[idx] = svlen
        ev_support[idx] = event.pe + event.supp + event.sc
        ev_svt[idx] = _c_svt_code(event.svtype)
        if tb < 0:
            continue  # chrB not present in bam; leave gt counts unset
        same = ta == tb
        map_records_by_tid[ta].append((max(0, start - 1), min(ref_lengths[ta], end + 1), idx))
        if not same:
            map_records_by_tid[tb].append((max(0, posB - 1), min(ref_lengths[tb], posB + 1), idx))
        if same and end - start <= GT_MAX_SPAN_ROI:
            rois_by_tid[ta].append((max(0, start - pad), min(ref_lengths[ta], end + pad)))
        else:
            rois_by_tid[ta].append((max(0, posA - pad), min(ref_lengths[ta], posA + pad)))
            rois_by_tid[tb].append((max(0, posB - pad), min(ref_lengths[tb], posB + pad)))

    cdef np.ndarray[np.int32_t, ndim=1] alt_counts = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.int32_t, ndim=1] ref_counts = np.zeros(n_events, dtype=np.int32)
    cdef np.ndarray[np.float64_t, ndim=1] ll_ref = np.zeros(n_events, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] ll_het = np.zeros(n_events, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] ll_hom = np.zeros(n_events, dtype=np.float64)
    cdef int* ev_ta_p = <int*>&ev_ta[0]
    cdef int* ev_tb_p = <int*>&ev_tb[0]
    cdef int* ev_posA_p = <int*>&ev_posA[0]
    cdef int* ev_posB_p = <int*>&ev_posB[0]
    cdef int* ev_start_p = <int*>&ev_start[0]
    cdef int* ev_end_p = <int*>&ev_end[0]
    cdef int* ev_span_p = <int*>&ev_span[0]
    cdef int* ev_svlen_p = <int*>&ev_svlen[0]
    cdef int* ev_support_p = <int*>&ev_support[0]
    cdef uint8_t* ev_svt_p = <uint8_t*>&ev_svt[0]
    cdef int* alt_counts_p = <int*>&alt_counts[0]
    cdef int* ref_counts_p = <int*>&ref_counts[0]
    cdef double* ll_ref_p = &ll_ref[0]
    cdef double* ll_het_p = &ll_het[0]
    cdef double* ll_hom_p = &ll_hom[0]
    cdef ankerl_set[uint64_t] alt_seen
    cdef ankerl_set[uint64_t] ref_seen
    cdef int clip_lo, clip_hi
    # Raw per-tid map records are shipped to scan jobs (in workers when W > 1);
    # each job builds its own si::IntervalMap in this insertion order, which is
    # exactly the order the Python wrapper used before, so candidate order and
    # downstream tie-breaking match earlier versions.
    ev_arrays = (ev_ta, ev_tb, ev_posA, ev_posB, ev_start, ev_end, ev_span, ev_svlen, ev_support, ev_svt)
    params = (mapq_threshold, min_sv_size, max_dist, insert_median, insert_stdev, clip_l, bool(paired_end), cslop)
    tid_order = list(rois_by_tid.keys())
    records = {tid: list(recs) for tid, recs in map_records_by_tid.items()}

    # breakpoint haplotypes for the junction scorer; disabled without a
    # reference (streamed input etc.) - votes then keep their flat weights
    ev_jx = None
    if paired_end and ref_genome is not None:
        try:
            ev_jx = _gt_build_junctions(events, ev_info, ref_genome, references)
        except Exception:
            logging.warning("Junction scorer disabled: reference fetch failed", exc_info=True)
            ev_jx = None

    # Shard the merged ROIs of each tid into W contiguous chunks. Chunks are
    # scanned independently (in workers when W > 1) and votes are replayed in
    # (tid, chunk) order, which reproduces single-pass scan order exactly.
    cdef int W = 1
    cdef int n_v = 0
    if procs > 1 and bam_path and n_events >= 2000:
        W = procs if procs < 8 else 8
    merged_by_tid = {}
    for tid in tid_order:
        intervals = sorted(rois_by_tid[tid])
        merged = []
        roi_start, roi_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= roi_end:
                roi_end = max(roi_end, e)
            else:
                merged.append((roi_start, roi_end))
                roi_start, roi_end = s, e
        merged.append((roi_start, roi_end))
        merged_by_tid[tid] = merged
    worker_units = []
    for w in range(W):
        units = {}
        for tid in tid_order:
            m = merged_by_tid[tid]
            s = len(m) * w // W
            e = len(m) * (w + 1) // W
            if e > s:
                units[tid] = m[s:e]
        worker_units.append(units)

    if W == 1:
        results = [_gt_scan_job(bam, tid_order, records, ev_arrays, params, worker_units[0], ev_jx)]
    else:
        payloads = [dict(bam_path=bam_path, tid_order=tid_order, records=records,
                         ev_arrays=ev_arrays, params=params, units=worker_units[w],
                         ev_jx=ev_jx) for w in range(W)]
        with multiprocessing.Pool(W) as pool:
            results = pool.map(_gt_worker, payloads)

    reads_seen = 0
    total_rois = 0
    n_alt_votes = 0
    n_ref_votes = 0
    n_jx_scored = 0
    n_jx_changed = 0
    for votes_j, rs_j, nr_j, js_j, jc_j in results:
        reads_seen += rs_j
        total_rois += nr_j
        n_jx_scored += js_j
        n_jx_changed += jc_j
    cdef np.ndarray[np.int32_t, ndim=1] v0
    cdef np.ndarray[np.uint64_t, ndim=1] v1
    cdef np.ndarray[np.uint8_t, ndim=1] v2
    cdef np.ndarray[np.float64_t, ndim=1] v3, v4
    for tid in tid_order:
        for w in range(W):
            vd = results[w][0].get(tid)
            if vd is None:
                continue
            v0 = vd[0]; v1 = vd[1]; v2 = vd[2]; v3 = vd[3]; v4 = vd[4]
            n_v = v0.shape[0]
            if n_v == 0:
                continue
            _gt_apply_votes(<int*>&v0[0], <uint64_t*>&v1[0], <uint8_t*>&v2[0],
                            &v3[0], &v4[0], n_v,
                            alt_counts_p, ref_counts_p, ll_ref_p, ll_het_p, ll_hom_p,
                            alt_seen, ref_seen, &c_alt, &c_ref)
            n_alt_votes += c_alt
            n_ref_votes += c_ref

    n_counted = 0
    for idx, event in enumerate(events):
        if ev_info[idx] is None or ev_info[idx]["tb"] < 0:
            continue
        if alt_counts[idx] == 0:
            continue  # no alt evidence found; leave unset so the default genotyping path is used
        event.gt_alt = alt_counts[idx]
        event.gt_ref = ref_counts[idx]
        event.gt_ll_ref = ll_ref[idx]
        event.gt_ll_het = ll_het[idx]
        event.gt_ll_hom = ll_hom[idx]
        n_counted += 1

    logging.info("Fragment genotyping: {} events ({} alt / {} ref votes)".format(n_events, n_alt_votes, n_ref_votes))
    logging.debug("Fragment genotyping detail: {:.2f}s, {} counted, {} ROIs, {} reads, junction scorer {} scored/{} replaced".format(
        time.time() - t0, n_counted, total_rois, reads_seen, n_jx_scored, n_jx_changed))
    return events


cdef int _c_classify(int f_tid, uint64_t qhash,
                     const int* cand, int n_cand,
                     const int* ev_ta, const int* ev_tb,
                     const int* ev_posA, const int* ev_posB,
                     const int* ev_start, const int* ev_end, const int* ev_span,
                     const int* ev_svlen, const uint8_t* ev_svt, const int* ev_support,
                     int f_start, int f_aend, int f_ext_lo, int f_ext_hi, int f_tlen,
                     bint f_proper, bint f_paired, int f_clip_lo, int f_clip_hi, int f_mtid, int f_mpos,
                     const int* op_v, const int* ln_v, const int* opos_v, int n_ops,
                     const int* c1s, const int* p1s, const int* c2s, const int* p2s, const int* qgs, int n_splits,
                     ankerl_set[uint64_t]& alt_seen, ankerl_set[uint64_t]& ref_seen,
                     bint paired_end, int max_dist, int insert_median, int insert_stdev,
                     int min_sv_size, int clip_l, int cslop,
                     int* out_ev, double* out_a, double* out_r) noexcept nogil:
    """Soft-likelihood assignment for a fragment: the best-matching candidate event
    receives an alt-weight contribution, otherwise the best-covered candidate receives
    a ref-weight contribution. A qname never contributes twice to the same event nor
    both alt and ref for it. Returns 1 (alt vote), 2 (ref vote) or 0; the winning
    event index and (a, r) weights are returned via out_ev/out_a/out_r so the caller
    can record the vote (see _gt_apply_votes for how votes become counts/LLs)."""
    cdef int k, idx, best_alt = -1, best_support = -1, cov_best_idx = -1, cov_best_sup = -1
    cdef double best_a = -1.0, best_r = 0.0, a, r, cov_a = 0.0, cov_r = 0.0
    cdef double oa = 0.0, orr = 0.0
    for k in range(n_cand):
        idx = cand[k]
        if ev_ta[idx] < 0:
            continue
        if f_tid != ev_ta[idx] and f_tid != ev_tb[idx]:
            continue
        _c_alt_weight(ev_ta[idx], ev_tb[idx], ev_posA[idx], ev_posB[idx], ev_start[idx], ev_end[idx],
                      ev_span[idx], ev_svlen[idx], <int>ev_svt[idx],
                      f_tid, f_start, f_aend, f_ext_lo, f_ext_hi, f_tlen, f_proper, f_paired,
                      f_clip_lo, f_clip_hi, f_mtid, f_mpos,
                      op_v, ln_v, opos_v, n_ops, c1s, p1s, c2s, p2s, qgs, n_splits,
                      paired_end, max_dist, insert_median, insert_stdev, min_sv_size, clip_l, cslop,
                      &oa, &orr)
        a = oa; r = orr
        if a >= 0 and (a > best_a or (a == best_a and ev_support[idx] > best_support)):
            best_alt = idx; best_a = a; best_r = r; best_support = ev_support[idx]
        if best_alt < 0:
            _c_ref_weight(ev_ta[idx], ev_tb[idx], ev_posA[idx], ev_posB[idx], ev_start[idx], ev_end[idx],
                          ev_span[idx], ev_svlen[idx], <int>ev_svt[idx],
                          f_start, f_aend, f_ext_lo, f_ext_hi, f_tlen, f_proper,
                          f_mtid, f_tid,
                          op_v, ln_v, opos_v, n_ops, c1s, p1s, c2s, p2s, qgs, n_splits,
                          paired_end, max_dist, insert_median, insert_stdev, min_sv_size, clip_l,
                          &oa, &orr)
            a = oa; r = orr
            if r >= 0:
                if cov_best_idx < 0 or ev_support[idx] > cov_best_sup:
                    cov_best_idx = idx; cov_a = a; cov_r = r; cov_best_sup = ev_support[idx]
    if best_alt >= 0:
        if ref_seen.find(support_key(best_alt, qhash)) == ref_seen.end() and \
                add_support_hash(alt_seen, best_alt, qhash):
            out_ev[0] = best_alt; out_a[0] = best_a; out_r[0] = best_r
            return 1
        return 0
    if cov_best_idx >= 0:
        # contribute to a single covered event: most independent support wins (stable choice)
        if alt_seen.find(support_key(cov_best_idx, qhash)) == alt_seen.end() and \
                add_support_hash(ref_seen, cov_best_idx, qhash):
            out_ev[0] = cov_best_idx; out_a[0] = cov_a; out_r[0] = cov_r
            return 2
    return 0


cdef void _gt_apply_votes(const int* v_ev, const uint64_t* v_qh, const uint8_t* v_kind,
                          const double* v_a, const double* v_r, int n,
                          int* alt_counts, int* ref_counts,
                          double* ll_ref, double* ll_het, double* ll_hom,
                          ankerl_set[uint64_t]& alt_seen, ankerl_set[uint64_t]& ref_seen,
                          int* out_alt, int* out_ref) noexcept nogil:
    # Replay votes in scan order through the global seen-sets: a vote is applied
    # only if this qname has not already voted (alt or ref) for the event. Given
    # votes in true scan order this reproduces a single-pass scan exactly.
    cdef int k, ev
    cdef uint64_t qh
    out_alt[0] = 0
    out_ref[0] = 0
    for k in range(n):
        ev = v_ev[k]
        qh = v_qh[k]
        if v_kind[k] == 1:
            if ref_seen.find(support_key(ev, qh)) == ref_seen.end() and \
                    add_support_hash(alt_seen, ev, qh):
                alt_counts[ev] += 1
                ll_ref[ev] += c_log(v_a[k] * 0.005 + v_r[k] * 0.995)
                ll_het[ev] += c_log(v_a[k] * 0.5 + v_r[k] * 0.5)
                ll_hom[ev] += c_log(v_a[k] * 0.99 + v_r[k] * 0.01)
                out_alt[0] += 1
        else:
            if alt_seen.find(support_key(ev, qh)) == alt_seen.end() and \
                    add_support_hash(ref_seen, ev, qh):
                ref_counts[ev] += 1
                ll_ref[ev] += c_log(v_a[k] * 0.005 + v_r[k] * 0.995)
                ll_het[ev] += c_log(v_a[k] * 0.5 + v_r[k] * 0.5)
                ll_hom[ev] += c_log(v_a[k] * 0.99 + v_r[k] * 0.01)
                out_ref[0] += 1


# --- Junction realignment scorer -------------------------------------------------
# Clean-room re-implementation of the manta-style junction lhood comparison
# (manta is PolyForm Strict; only the mathematical idea and parameter values are
# used). For DEL/INS events with known inserted sequence we build the two
# breakpoint haplotype strings: J_ref = reference across the junction, J_alt =
# reference with the allele applied (deleted span removed / inserted sequence
# added). A read's junction-spanning segment (collected from its CIGAR in
# J_ref coordinates, reverse-complemented for minus strand) is compared
# ungapped against both haplotypes over a small offset scan; per-base log
# probabilities come from base qualities blended with a SNP-aware error
# (theta=1e-3, manta's value). The scan replaces the flat readwin-ref
# (0.02, 0.95) and INS clip (0.75, 0.10) vote weights with per-fragment
# lhood-derived (a, r): absorbed repeat-insertion reads then score alt instead
# of ref, and clips matching the inserted sequence are upgraded while clips
# that match neither junction are neutralised. Manta's parameters: >=16 bp
# flank per side, mismatch density < 0.25, mismap floors 1e-5 alt / 1e-6 ref.
cdef int GT_JX_FLANK = 150
cdef int GT_JX_MAX_SVLEN = 1000
cdef int GT_JX_MIN_SIDE = 16
cdef int GT_JX_SCAN = 24
cdef double GT_JX_UP = -1.3862943611198906  # log(0.25): per-base cost of an unaligned base
cdef double GT_JX_FLOOR_A = 1e-5
cdef double GT_JX_FLOOR_R = 1e-6

cdef double GT_JX_LM[128]  # log P(observed base matches | BQ)
cdef double GT_JX_LF[128]  # log P(mismatch | BQ)
cdef unsigned char GT_JX_BCODE[256]
cdef unsigned char GT_JX_BCOMP[256]
_GT_BCODE_TAB = b""  # bytes translation table (set by _init_gt_jx_tables)
_GT_ACODE_TAB = b""   # code -> ASCII (set by _init_gt_jx_tables)


def _init_gt_jx_tables():
    from math import log
    global _GT_BCODE_TAB, _GT_ACODE_TAB
    cdef int q, i
    cdef double p
    for q in range(128):
        p = 1e-3 + 0.999 * 10.0 ** (-q / 10.0)
        if p > 1.0 - 1e-15:
            p = 1.0 - 1e-15
        GT_JX_LM[q] = log(1.0 - p)
        GT_JX_LF[q] = log(p)
    for i in range(256):
        GT_JX_BCODE[i] = 15
        GT_JX_BCOMP[i] = 15
    for ch, code in ((65, 0), (67, 1), (71, 2), (84, 3)):  # A C G T
        GT_JX_BCODE[ch] = code
        GT_JX_BCODE[ch + 32] = code
    GT_JX_BCOMP[0] = 3
    GT_JX_BCOMP[1] = 2
    GT_JX_BCOMP[2] = 1
    GT_JX_BCOMP[3] = 0
    _GT_BCODE_TAB = bytes(bytearray(GT_JX_BCODE[i] for i in range(256)))
    t = bytearray(256)
    for i in range(256):
        t[i] = 78  # N
    t[0] = 65
    t[1] = 67
    t[2] = 71
    t[3] = 84
    _GT_ACODE_TAB = bytes(t)


_init_gt_jx_tables()


def _gt_extract_ins_seq(event, int svlen):
    # inserted sequence for the alt haplotype: prefer the assembled variant
    # sequence, else the lowercase (unaligned) middle of the contig between
    # its two aligned blocks
    vs = event.variant_seq
    if isinstance(vs, str) and not vs.startswith("<"):
        if len(vs) >= max(15, svlen // 2):
            return vs.upper()
    contig = event.contig
    if isinstance(contig, str) and event.contig_lc > 0 and event.contig_rc > 0:
        mid = contig[event.contig_lc: len(contig) - event.contig_rc]
        if len(mid) >= max(15, svlen // 2) and \
                sum(1 for ch in mid if ch.islower()) * 2 >= len(mid):
            return mid.upper()
    return None


def _gt_build_junctions(events, ev_info, ref_genome, bam_refs):
    # one (kind, J_ref, J_alt, F, del_len, ins_len, posA, posB) tuple per
    # eligible event, else None. J strings are upper-case reference sequence;
    # coordinates are BAM 0-based (same space as event.posA/posB).
    ref_lens = {}
    for rname, rlen in zip(ref_genome.references, ref_genome.lengths):
        ref_lens[rname] = rlen
    by_chrom = {}
    for br in bam_refs:
        if br in ref_lens:
            by_chrom[br] = br
        elif br.startswith("chr") and br[3:] in ref_lens:
            by_chrom[br] = br[3:]
        elif "chr" + br in ref_lens:
            by_chrom[br] = "chr" + br
    cdef int F = GT_JX_FLANK
    # events per chromosome with junction windows, for the complex-locus
    # guard: the scorer assumes a single junction inside its window
    win_by_chrom = {}
    for idx, info in enumerate(ev_info):
        if info is None or info["tb"] < 0 or info["ta"] != info["tb"]:
            continue
        if info["svtype"] not in ("DEL", "INS"):
            continue
        posA, posB = info["posA"], info["posB"]
        lo = posA - F
        hi = posB + F
        win_by_chrom.setdefault(info["ta"], []).append((lo, hi, idx))
    for k in win_by_chrom:
        win_by_chrom[k].sort()
    jx = []
    n_del = n_ins = 0
    for idx, event in enumerate(events):
        info = ev_info[idx]
        if info is None or info["tb"] < 0 or info["ta"] != info["tb"]:
            jx.append(None)
            continue
        svt = info["svtype"]
        if svt not in ("DEL", "INS"):
            jx.append(None)
            continue
        posA, posB = info["posA"], info["posB"]
        lo = posA - F
        hi = posB + F
        wins = win_by_chrom.get(info["ta"], [])
        n_others = 0
        for wlo, whi, widx in wins:
            if wlo > hi:
                break
            if widx != idx and whi >= lo:
                n_others += 1
                if n_others >= 2:
                    break
        if n_others >= 2:
            jx.append(None)
            continue
        chrom = by_chrom.get(bam_refs[info["ta"]])
        if chrom is None:
            jx.append(None)
            continue
        if svt == "DEL":
            if posB - posA < 1 or posB - posA > GT_JX_MAX_SVLEN:
                jx.append(None)
                continue
            lo, hi = posA - F, posB + F
            if lo < 0 or hi > ref_lens[chrom]:
                jx.append(None)
                continue
            J_ref = ref_genome.fetch(chrom, lo, hi).upper()
            J_alt = ref_genome.fetch(chrom, lo, posA).upper() + ref_genome.fetch(chrom, posB, hi).upper()
            jx.append((0, J_ref, J_alt, F, posB - posA, 0, posA, posB))
            n_del += 1
        else:
            if i_abs(posB - posA) > 8:
                jx.append(None)
                continue
            anchor = (posA + posB) // 2
            svlen = info["svlen"]
            ins = _gt_extract_ins_seq(event, svlen)
            lo, hi = anchor - F, anchor + F
            if ins is None or lo < 0 or hi > ref_lens[chrom]:
                jx.append(None)
                continue
            J_ref = ref_genome.fetch(chrom, lo, hi).upper()
            J_alt = ref_genome.fetch(chrom, lo, anchor).upper() + ins + ref_genome.fetch(chrom, anchor, hi).upper()
            jx.append((1, J_ref, J_alt, F, 0, len(ins), anchor, anchor))
            n_ins += 1
    logging.debug("Junction scorer: {} DEL + {} INS junctions built".format(n_del, n_ins))
    return jx


cdef void _c_scan_best(const uint8_t* sb, const uint8_t* sq, const int* sc, const int n,
                       const uint8_t* jb, const int m, const int F,
                       const int del_len, const int ins_len,
                       const int o_lo, const int o_hi,
                       const bint check_flank, const int jpos,
                       double* out_l, int* out_ok) noexcept nogil:
    # Best ungapped BQ-weighted log-lhood of the read junction-string S
    # against junction haplotype J over offsets o in [o_lo, o_hi]. S coords
    # are J_ref-space; mapping to this allele's space: c < F -> c;
    # F <= c < F+del_len -> unaligned (deleted span vs alt haplotype);
    # else c + ins_len - del_len. Bases that map outside [0, m) or are
    # skipped cost GT_JX_UP each. Junction bases that lie inside the read's
    # aligned span but have no read base covering them (the inserted region
    # of J_alt unseen by a ref read, or the deleted span of J_ref unseen by
    # an alt read) cost GT_JX_UP each - without this the two haplotypes tie
    # for any read that does not physically span the novel junction.
    # Offsets failing the flank rule (>= GT_JX_MIN_SIDE aligned bases per
    # side of the junction point jpos), the overlap floor (2*GT_JX_MIN_SIDE)
    # or the mismatch-density filter (> 0.25) are rejected.
    cdef int o, i, c, p, ov, mm, left, right, bq, p_min, p_max
    cdef double L, best = -1e18
    cdef bint any_ok = False
    for o in range(o_lo, o_hi + 1):
        L = 0.0
        ov = 0
        mm = 0
        left = 0
        right = 0
        p_min = 2147483647
        p_max = -1
        for i in range(n):
            c = sc[i]
            if c >= F:
                if del_len and c < F + del_len:
                    L += GT_JX_UP
                    continue
                p = c + ins_len - del_len
            else:
                p = c
            p += o
            if p < 0 or p >= m:
                L += GT_JX_UP
                continue
            bq = sq[i]
            if bq > 93:
                bq = 93
            if sb[i] == jb[p] and sb[i] != 15:
                L += GT_JX_LM[bq]
            else:
                L += GT_JX_LF[bq]
                mm += 1
            ov += 1
            if p < p_min:
                p_min = p
            if p > p_max:
                p_max = p
            if check_flank:
                if p < jpos:
                    left += 1
                else:
                    right += 1
        if ov == 0:
            continue
        L += <double>(p_max - p_min + 1 - ov) * GT_JX_UP
        if check_flank and (left < GT_JX_MIN_SIDE or right < GT_JX_MIN_SIDE):
            continue
        if ov < 2 * GT_JX_MIN_SIDE:
            continue
        if mm * 4 > ov:
            continue
        if L > best:
            best = L
            any_ok = True
    out_l[0] = best
    out_ok[0] = 1 if any_ok else 0


cdef inline void _c_jx_weights(double la, double lr, double* out_a, double* out_r) noexcept nogil:
    # lhoods -> (a, r) vote weights: normalise with manta's mismap floors, then
    # clamp to keep a fragment of garbage evidence from being decisive
    cdef double M, za, zr, pa
    M = la if la > lr else lr
    za = c_exp(la - M)
    zr = c_exp(lr - M)
    if za < GT_JX_FLOOR_A:
        za = GT_JX_FLOOR_A
    if zr < GT_JX_FLOOR_R:
        zr = GT_JX_FLOOR_R
    pa = za / (za + zr)
    if pa < 0.02:
        pa = 0.02
    elif pa > 0.98:
        pa = 0.98
    out_a[0] = pa
    out_r[0] = 1.0 - pa


def _gt_gapped_lhood(s_ascii, const uint8_t[:] sqv, str j_seq):
    # edlib HW realignment of the read junction-string against one haplotype;
    # BQ-weighted lhood from the extended cigar. Gaps and target-only bases
    # cost GT_JX_UP per base, matching the ungapped span-gap model. Returns
    # (lhood, ok).
    from dysgu.edlib import edlib
    cdef int n = len(s_ascii)
    cdef int k = n // 3
    if k < 20:
        k = 20
    r1 = edlib.align(s_ascii, j_seq, mode="HW", task="path", k=k)
    if r1["editDistance"] < 0:
        return -1e18, False
    cdef double L = 0.0
    cdef int qi = 0, num = 0
    cdef char op
    cdef int bq
    cigar = r1["cigar"]
    cdef bytes cb = cigar.encode() if isinstance(cigar, str) else cigar
    cdef const char* cp = cb
    cdef int i = 0
    cdef int cl = len(cb)
    while i < cl:
        if cp[i] >= 48 and cp[i] <= 57:
            num = num * 10 + (cp[i] - 48)
            i += 1
            continue
        op = cp[i]
        i += 1
        if op == 61:  # '='
            for _ in range(num):
                bq = sqv[qi]
                if bq > 93:
                    bq = 93
                L += GT_JX_LM[bq]
                qi += 1
        elif op == 88:  # 'X'
            for _ in range(num):
                bq = sqv[qi]
                if bq > 93:
                    bq = 93
                L += GT_JX_LF[bq]
                qi += 1
        elif op == 73:  # 'I': query bases absent from target
            for _ in range(num):
                L += GT_JX_UP
                qi += 1
        else:  # 'D': target bases unobserved by the query
            L += <double>num * GT_JX_UP
        num = 0
    return L, True


def _gt_junction_span(aln, tuple jx, int read_start):
    # Score a spanning read against both breakpoint haplotypes; returns
    # (a, r) weights or (-1.0, -1.0) when the read is not junction-informative.
    cdef int kind = jx[0]
    cdef int F = jx[3], del_len = jx[4], ins_len = jx[5]
    cdef int posA = jx[6], posB = jx[7]
    seq = aln.query_sequence
    if seq is None:
        return -1.0, -1.0
    cdef int wlo, whi
    if kind == 0:
        wlo = posA - F
        whi = posB + F
    else:
        wlo = posA - F
        whi = posA + F
    if whi - wlo > 4096:
        return -1.0, -1.0
    s8 = seq.encode()
    q = aln.query_qualities
    cdef bint have_q = q is not None
    sb = bytearray()
    sq = bytearray()
    coords = []
    cdef int qpos = 0, rpos = read_start, op, ln, k, c
    for op, ln in aln.cigartuples:
        if op == 0 or op == 7 or op == 8:  # M = X
            for k in range(ln):
                if wlo <= rpos < whi:
                    sb.append(GT_JX_BCODE[s8[qpos]])
                    sq.append(q[qpos] if have_q else 30)
                    coords.append(rpos - wlo)
                qpos += 1
                rpos += 1
        elif op == 1:  # I: inserted bases extend the anchor's coord
            c = rpos - 1
            for k in range(ln):
                if wlo <= c < whi:
                    sb.append(GT_JX_BCODE[s8[qpos]])
                    sq.append(q[qpos] if have_q else 30)
                    coords.append(c - wlo + 1 + k)
                qpos += 1
        elif op == 2 or op == 3:  # D, N
            rpos += ln
        else:  # S, H, P
            qpos += ln
    cdef int n = len(sb)
    if n < 2 * GT_JX_MIN_SIDE + 2:
        return -1.0, -1.0
    # query_sequence is +strand-oriented for both strands (BAM stores SEQ
    # reverse-complemented for reverse alignments), so the cigar walk above
    # already yields bases in reference-forward order - no RC needed
    cdef const uint8_t[:] sbv = sb
    cdef const uint8_t[:] sqv = sq
    cdef np.ndarray[np.int32_t, ndim=1] scv = np.array(coords, dtype=np.int32)
    J_ref_b = jx[1].encode()
    J_alt_b = jx[2].encode()
    cdef const uint8_t[:] jrv = J_ref_b.translate(_GT_BCODE_TAB)
    cdef const uint8_t[:] jav = J_alt_b.translate(_GT_BCODE_TAB)
    cdef double la = 0.0, lr = 0.0
    cdef int oka = 0, okr = 0
    cdef int m_ref = len(jrv), m_alt = len(jav)
    cdef int jpos = F + ins_len
    with nogil:
        _c_scan_best(&sbv[0], &sqv[0], <const int*>&scv[0], n, &jav[0], m_alt, F, del_len, ins_len,
                     -GT_JX_SCAN, GT_JX_SCAN, True, jpos, &la, &oka)
        _c_scan_best(&sbv[0], &sqv[0], <const int*>&scv[0], n, &jrv[0], m_ref, F, 0, 0,
                     -GT_JX_SCAN, GT_JX_SCAN, False, jpos, &lr, &okr)
    if not oka or not okr:
        return -1.0, -1.0
    cdef double na = 0.0, nr = 0.0
    _c_jx_weights(la, lr, &na, &nr)
    if kind == 1 and 0.05 < na < 0.95:
        # ungapped placement is ambiguous (repeat-absorbed read): fall back
        # to a gapped HW realignment against both haplotypes, which can place
        # the absorbed bases into the inserted region of J_alt. DEL needs no
        # gapped fallback - its signal is the ungapped span-gap penalty
        s_ascii = bytes(sb).translate(_GT_ACODE_TAB)
        la_g, okg_a = _gt_gapped_lhood(s_ascii, sqv, jx[2])
        lr_g, okg_r = _gt_gapped_lhood(s_ascii, sqv, jx[1])
        if okg_a and okg_r:
            _c_jx_weights(la_g, lr_g, &na, &nr)
    return na, nr


def _gt_junction_clip(aln, tuple jx):
    # Score a soft-clipped read's clip against the inserted sequence and the
    # reference continuation on the clipped side; returns (a, r) or (-1, -1).
    cdef int F = jx[3], ins_len = jx[5]
    if ins_len < 8:
        return -1.0, -1.0
    seq = aln.query_sequence
    if seq is None:
        return -1.0, -1.0
    s8 = seq.encode()
    q = aln.query_qualities
    cdef bint have_q = q is not None
    ct = aln.cigartuples
    cdef list clips = []  # (codes bytearray, quals bytearray, is_leading)
    if ct[0][0] == 4 and ct[0][1] >= 8:  # leading S
        ln = ct[0][1]
        clips.append((bytearray(s8[:ln]), bytearray(q[:ln]) if have_q else bytearray([30] * ln), True))
    if ct[-1][0] == 4 and ct[-1][1] >= 8:  # trailing S
        ln = ct[-1][1]
        clips.append((bytearray(s8[len(s8) - ln:]), bytearray(q[len(s8) - ln:]) if have_q else bytearray([30] * ln), False))
    if not clips:
        return -1.0, -1.0
    J_alt_b = jx[2].encode()
    cdef const uint8_t[:] jav = J_alt_b.translate(_GT_BCODE_TAB)
    # candidates: inserted sequence (either strand) and the reference
    # continuation on the clipped side (left of breakpoint for a 5' clip,
    # right of it for a 3' clip)
    ins_b = bytes(jav[F: F + ins_len])
    ins_rc = bytes([GT_JX_BCOMP[b] for b in reversed(ins_b)])
    left_b = bytes(jav[:F])
    right_b = bytes(jav[F + ins_len:])
    cdef const uint8_t[:] ins_v = ins_b
    cdef const uint8_t[:] insr_v = ins_rc
    cdef const uint8_t[:] left_v = left_b
    cdef const uint8_t[:] right_v = right_b
    cdef int ins_n = len(ins_b)
    cdef double la = -1e18, lr = -1e18, l1, l2
    cdef int ok = 0
    cdef int cn, o_lo
    cdef const uint8_t[:] cv
    cdef const uint8_t[:] qv
    cdef np.ndarray[np.int32_t, ndim=1] dummy
    cdef const uint8_t* ref_p
    cdef int ref_n
    for cb, cq, leading in clips:
        cn = len(cb)
        o_lo = -cn + 1
        cv = cb
        qv = cq
        dummy = np.arange(cn, dtype=np.int32)
        if leading:
            ref_p = &left_v[0]
            ref_n = len(left_b)
        else:
            ref_p = &right_v[0]
            ref_n = len(right_b)
        with nogil:
            _c_scan_best(&cv[0], &qv[0], <const int*>&dummy[0], cn, &ins_v[0], ins_n, F, 0, 0,
                         o_lo, ins_n - 1, False, 0, &l1, &ok)
            if ok and l1 > la:
                la = l1
            _c_scan_best(&cv[0], &qv[0], <const int*>&dummy[0], cn, &insr_v[0], ins_n, F, 0, 0,
                         o_lo, ins_n - 1, False, 0, &l2, &ok)
            if ok and l2 > la:
                la = l2
            _c_scan_best(&cv[0], &qv[0], <const int*>&dummy[0], cn, ref_p, ref_n, F, 0, 0,
                         o_lo, ref_n - 1, False, 0, &l1, &ok)
            if ok and l1 > lr:
                lr = l1
    if la <= -1e17 or lr <= -1e17:
        return -1.0, -1.0
    cdef double na = 0.0, nr = 0.0
    _c_jx_weights(la, lr, &na, &nr)
    if 0.25 < na < 0.9:
        # clip matches the inserted sequence and the reference continuation
        # about equally (tandem repeat): the sequence content is genuinely
        # uninformative - keep the original cigar-level clip vote, which
        # carries the breakpoint-position signal
        return -1.0, -1.0
    return na, nr


def _gt_scan_job(object bam, list tid_order, dict records, tuple ev_arrays, tuple params, dict units, list ev_jx=None):
    """Scan this job's share of merged ROIs per tid and classify fetched reads against
    the global per-tid map records. Votes (event, qhash, kind, a, r) are returned per
    tid in scan order together with (reads_seen, n_rois); the caller replays the votes
    through the global seen-sets, so results are identical to a single-pass scan no
    matter how the ROIs are sharded across jobs."""
    from dysgu.graph import AlignmentsSA
    cdef int mapq_t, min_sv_t, max_dist_t, ins_med, ins_sd, clip_lt, cslop, tid_c
    cdef bint paired_end
    mapq_t, min_sv_t, max_dist_t, ins_med, ins_sd, clip_lt, paired_end, cslop = params
    cdef np.ndarray[np.int32_t, ndim=1] ev_ta = ev_arrays[0]
    cdef np.ndarray[np.int32_t, ndim=1] ev_tb = ev_arrays[1]
    cdef np.ndarray[np.int32_t, ndim=1] ev_posA = ev_arrays[2]
    cdef np.ndarray[np.int32_t, ndim=1] ev_posB = ev_arrays[3]
    cdef np.ndarray[np.int32_t, ndim=1] ev_start = ev_arrays[4]
    cdef np.ndarray[np.int32_t, ndim=1] ev_end = ev_arrays[5]
    cdef np.ndarray[np.int32_t, ndim=1] ev_span = ev_arrays[6]
    cdef np.ndarray[np.int32_t, ndim=1] ev_svlen = ev_arrays[7]
    cdef np.ndarray[np.int32_t, ndim=1] ev_support = ev_arrays[8]
    cdef np.ndarray[np.uint8_t, ndim=1] ev_svt = ev_arrays[9]
    cdef int* ev_ta_p = <int*>&ev_ta[0]
    cdef int* ev_tb_p = <int*>&ev_tb[0]
    cdef int* ev_posA_p = <int*>&ev_posA[0]
    cdef int* ev_posB_p = <int*>&ev_posB[0]
    cdef int* ev_start_p = <int*>&ev_start[0]
    cdef int* ev_end_p = <int*>&ev_end[0]
    cdef int* ev_span_p = <int*>&ev_span[0]
    cdef int* ev_svlen_p = <int*>&ev_svlen[0]
    cdef int* ev_support_p = <int*>&ev_support[0]
    cdef uint8_t* ev_svt_p = <uint8_t*>&ev_svt[0]
    cdef int start, aend, pos, op, ln, mtid, mpos, tlen, ext_lo, ext_hi, lo, i
    cdef int cigar_value, n_cand, vote, clip_lo, clip_hi, out_ev
    cdef int n_jx_scored = 0, n_jx_changed = 0
    cdef long long isize
    cdef bint proper, paired
    cdef uint32_t cigar_l
    cdef uint32_t *cigar_p
    cdef uint64_t qhash
    cdef double out_a, out_r
    cdef AlignedSegment aln
    cdef cpp_vector[int] v_op, v_oln, v_opos
    cdef cpp_vector[int] v_found
    cdef cpp_vector[int] v_c1, v_p1, v_c2, v_p2, v_qg
    cdef cpp_vector[int] v_vev
    cdef cpp_vector[uint8_t] v_vkind
    cdef cpp_vector[uint64_t] v_vqh
    cdef cpp_vector[double] v_va, v_vr
    cdef int[::1] mv_ev
    cdef uint64_t[::1] mv_qh
    cdef uint8_t[::1] mv_kd
    cdef double[::1] mv_a, mv_r
    cdef ankerl_set[uint64_t] alt_seen
    cdef ankerl_set[uint64_t] ref_seen
    cdef SiIntMap simap
    reads_seen = 0
    n_rois = 0
    votes = {}
    for tid in tid_order:
        roi_list = units.get(tid)
        if not roi_list:
            continue
        rec_list = records.get(tid)
        if not rec_list:
            continue
        tid_c = tid
        # rebuild the interval index in the original map-record insertion order;
        # si::IntervalMap sorts internally exactly as the Python wrapper did, so
        # candidate order (and downstream tie-breaking) matches earlier versions
        simap.clear()
        for ms0, me0, idx0 in rec_list:
            simap.add(ms0, me0, idx0)
        simap.build()
        v_vev.clear(); v_vqh.clear(); v_vkind.clear(); v_va.clear(); v_vr.clear()
        for rs, re_ in roi_list:
            n_rois += 1
            for aln in bam.fetch(bam.references[tid], rs, re_):
                if aln._delegate.core.flag & 3844:
                    continue
                if aln._delegate.core.qual < mapq_t:
                    continue
                cigar_l = aln._delegate.core.n_cigar
                if cigar_l == 0:
                    continue
                start = aln._delegate.core.pos
                cigar_p = bam_get_cigar(aln._delegate)
                v_op.clear(); v_oln.clear(); v_opos.clear()
                pos = start
                aend = start
                for i in range(cigar_l):
                    cigar_value = <int>cigar_p[i]
                    op = cigar_value & BAM_CIGAR_MASK
                    ln = cigar_value >> BAM_CIGAR_SHIFT
                    if op == CINS or op == CDEL or op == CREF_SKIP:
                        v_op.push_back(op); v_oln.push_back(ln); v_opos.push_back(pos)
                    if op == CMATCH or op == CDEL or op == CREF_SKIP or op == CEQUAL or op == CDIFF:
                        pos += ln
                        aend += ln
                    elif op == CSOFT_CLIP or op == CHARD_CLIP:
                        v_op.push_back(op); v_oln.push_back(ln); v_opos.push_back(pos)
                reads_seen += 1
                mtid = aln._delegate.core.mtid
                mpos = aln._delegate.core.mpos
                proper = (aln._delegate.core.flag & 2) != 0
                isize = aln._delegate.core.isize
                tlen = <int>(isize if isize >= 0 else -isize) if isize else 0
                paired = (aln._delegate.core.flag & 1) != 0 and (aln._delegate.core.flag & 8) == 0 and mpos >= 0
                ext_lo = start
                ext_hi = aend
                if paired_end and proper and tlen > 0 and mtid == tid_c:
                    # exact pair span for concordant PE fragments
                    lo = min(start, mpos)
                    ext_lo = lo
                    ext_hi = lo + tlen
                elif paired and mtid == tid_c and mpos >= 0:
                    ext_lo = min(start, mpos)
                    ext_hi = max(aend, mpos)
                v_found.clear()
                simap.search_values(ext_lo, ext_hi + 1, v_found)
                n_cand = <int>v_found.size()
                if n_cand == 0:
                    continue
                clip_lo = 0
                clip_hi = 0
                clip_sizes_hard(aln, &clip_lo, &clip_hi)
                v_c1.clear(); v_p1.clear(); v_c2.clear(); v_p2.clear(); v_qg.clear()
                if aln.has_tag("SA"):
                    try:
                        sa = AlignmentsSA(aln, bam.gettid, paired_end)
                        sa.connect_alignments(aln, max_dist_t, mapq_t, 1)
                        for j in sa.join_result:
                            v_c1.push_back(j.chrom); v_p1.push_back(j.pos)
                            v_c2.push_back(j.chrom2); v_p2.push_back(j.pos2)
                            v_qg.push_back(j.query_end - j.query_pos)
                    except Exception:
                        v_c1.clear(); v_p1.clear(); v_c2.clear(); v_p2.clear(); v_qg.clear()
                qhash = aln_qhash(aln)
                with nogil:
                    vote = _c_classify(tid_c, qhash, v_found.data(), n_cand,
                                       ev_ta_p, ev_tb_p, ev_posA_p, ev_posB_p, ev_start_p, ev_end_p,
                                       ev_span_p, ev_svlen_p, ev_svt_p, ev_support_p,
                                       start, aend, ext_lo, ext_hi, tlen, proper, paired,
                                       clip_lo, clip_hi, mtid, mpos,
                                       v_op.data(), v_oln.data(), v_opos.data(), <int>v_op.size(),
                                       v_c1.data(), v_p1.data(), v_c2.data(), v_p2.data(), v_qg.data(), <int>v_qg.size(),
                                       alt_seen, ref_seen,
                                       paired_end, max_dist_t, ins_med, ins_sd, min_sv_t, clip_lt, cslop,
                                       &out_ev, &out_a, &out_r)
                if vote == 1:
                    if ev_jx is not None and out_a == 0.75 and out_r == 0.10:
                        # INS breakend clip vote: rescore the clip against the
                        # inserted sequence and the reference continuation
                        jx = ev_jx[out_ev]
                        if jx is not None:
                            n_jx_scored += 1
                            na, nr = _gt_junction_clip(aln, jx)
                            if nr >= 0:
                                out_a, out_r = na, nr
                                n_jx_changed += 1
                    v_vev.push_back(out_ev); v_vqh.push_back(qhash); v_vkind.push_back(1)
                    v_va.push_back(out_a); v_vr.push_back(out_r)
                elif vote == 2:
                    if ev_jx is not None and out_a == 0.02 and out_r == 0.95:
                        # read-level window cover vote: rescore the read's
                        # junction span against both haplotypes
                        jx = ev_jx[out_ev]
                        if jx is not None:
                            n_jx_scored += 1
                            na, nr = _gt_junction_span(aln, jx, start)
                            if nr >= 0:
                                out_a, out_r = na, nr
                                n_jx_changed += 1
                    v_vev.push_back(out_ev); v_vqh.push_back(qhash); v_vkind.push_back(2)
                    v_va.push_back(out_a); v_vr.push_back(out_r)
        if v_vev.empty():
            votes[tid] = (np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.uint64),
                          np.zeros(0, dtype=np.uint8), np.zeros(0, dtype=np.float64),
                          np.zeros(0, dtype=np.float64))
        else:
            mv_ev = <int[:v_vev.size()]>v_vev.data()
            mv_qh = <uint64_t[:v_vqh.size()]>v_vqh.data()
            mv_kd = <uint8_t[:v_vkind.size()]>v_vkind.data()
            mv_a = <double[:v_va.size()]>v_va.data()
            mv_r = <double[:v_vr.size()]>v_vr.data()
            # np.array copies: raw pointer-cast memoryviews export zero-copy views
            # whose backing vector storage is reused on the next tid
            votes[tid] = (np.array(mv_ev), np.array(mv_qh), np.array(mv_kd),
                          np.array(mv_a), np.array(mv_r))
    return votes, reads_seen, n_rois, n_jx_scored, n_jx_changed


def _gt_worker(dict payload):
    bam = pysam.AlignmentFile(payload["bam_path"], "rb")
    try:
        return _gt_scan_job(bam, payload["tid_order"], payload["records"],
                            payload["ev_arrays"], payload["params"], payload["units"],
                            payload.get("ev_jx"))
    finally:
        bam.close()

# https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
cdef int[32] tab32 = [
    0,  9,  1, 10, 13, 21,  2, 29,
    11, 14, 16, 18, 22, 25,  3, 30,
     8, 12, 20, 28, 15, 17, 24,  7,
    19, 27, 23,  6, 26,  5,  4, 31]

cdef int log2_32(uint32_t value):
    value |= value >> 1
    value |= value >> 2
    value |= value >> 4
    value |= value >> 8
    value |= value >> 16
    return tab32[<uint32_t>(value*0x07C4ACDD) >> 27]


cdef void window_rate(WindowRate *result, uint32_t cigar_l, uint32_t *cigar_p, int index, int window_size, bint reverse):
    cdef int n = 0
    cdef int matches = 0
    cdef int covered = 0
    cdef int i = index
    cdef int stop, opp, l
    cdef uint32_t cigar_value
    if reverse:
        stop = 0
    else:
        stop = cigar_l - 1
    while i != stop:
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & BAM_CIGAR_MASK
        l = <int> cigar_value >> BAM_CIGAR_SHIFT
        if opp == CMATCH or opp == CEQUAL:
            matches += l
            covered += l
        elif opp == CDEL:
            covered += l
            n += 1
        elif opp == CINS or opp == CSOFT_CLIP or opp == CHARD_CLIP or opp == CDIFF:
            n += 1
        if reverse:
            i -= 1
        else:
            i += 1
        if covered >= window_size:
            # if final match block is really long, the previous cigar opp can be missed if this is not done:
            matches = min(window_size, matches)
            break
    result.rate = 0 if not matches else n / (matches + n)
    result.index = i


cpdef filter_poorly_aligned_ends(spanning_alignments, float divergence=0.02):
    # This trys to check if the ends of reads are poorly aligned by looking at the rate of non-matching cigar events
    cdef int i, j, index_begin, index_end, cigar_index, stop_idx
    cdef int ws = 2000
    cdef int opp, length, cov_acc, mat_acc, n_acc, mat_acc_capped
    cdef int matches_total, n_total
    cdef float threshold, w_rate, mean, std
    cdef uint32_t cigar_value
    cdef uint32_t cigar_l
    cdef uint32_t *cigar_p
    cdef AlignedSegment alignment
    cdef list cov_ops, mat_ops, n_ops
    rates = []
    for item in spanning_alignments:
        alignment = item.align
        cigar_l = alignment._delegate.core.n_cigar
        if cigar_l == 0:
            continue
        cigar_p = bam_get_cigar(alignment._delegate)
        matches_total = 0
        n_total = 0
        for i in range(cigar_l):
            cigar_value = cigar_p[i]
            opp = <int> cigar_value & BAM_CIGAR_MASK
            length = <int> cigar_value >> BAM_CIGAR_SHIFT
            if opp == CMATCH or opp == CEQUAL:
                matches_total += length
            elif opp == CDEL or opp == CINS or opp == CSOFT_CLIP or opp == CHARD_CLIP or opp == CDIFF:
                n_total += 1
        if matches_total == 0:
            continue  # no aligned bases; rate undefined, exclude from statistics
        else:
            rates.append(n_total / float(matches_total + n_total))
    if rates:
        mean = np.mean(rates)
        std = np.std(rates)
    else:
        mean = 0
        std = 0

    multiplier = 3
    threshold = min(divergence, max(0.008, multiplier*mean, mean + (std*3)))
    spanning = []
    for item in spanning_alignments:
        alignment = item.align
        cigar_l = alignment._delegate.core.n_cigar
        if cigar_l == 0:
            continue
        cigar_p = bam_get_cigar(alignment._delegate)
        cigar_index = item.cigar_index

        cov_ops = [0] * <int> cigar_l
        mat_ops = [0] * <int> cigar_l
        n_ops = [0] * <int> cigar_l
        for i in range(cigar_l):
            cigar_value = cigar_p[i]
            opp = <int> cigar_value & BAM_CIGAR_MASK
            length = <int> cigar_value >> BAM_CIGAR_SHIFT
            if opp == CMATCH or opp == CEQUAL:
                cov_ops[i] = length
                mat_ops[i] = length
            elif opp == CDEL:
                cov_ops[i] = length
                n_ops[i] = 1
            elif opp == CINS or opp == CSOFT_CLIP or opp == CHARD_CLIP or opp == CDIFF:
                n_ops[i] = 1

        # Forward scan: find the first poorly-aligned prefix window.
        index_begin = 0
        j = 0
        cov_acc = 0
        mat_acc = 0
        n_acc = 0
        for i in range(<int> cigar_l):
            while j < <int> cigar_l and cov_acc < ws:
                cov_acc += cov_ops[j]
                mat_acc += mat_ops[j]
                n_acc += n_ops[j]
                j += 1
            if cov_acc >= ws:
                mat_acc_capped = mat_acc if mat_acc < ws else ws
                if mat_acc_capped == 0:
                    w_rate = 0.0
                else:
                    w_rate = n_acc / float(mat_acc_capped + n_acc)
                stop_idx = j
            else:
                if mat_acc == 0:
                    w_rate = 0.0
                else:
                    w_rate = n_acc / float(mat_acc + n_acc)
                stop_idx = <int> cigar_l
            if w_rate < threshold:
                break
            index_begin = stop_idx
            cov_acc -= cov_ops[i]
            mat_acc -= mat_ops[i]
            n_acc -= n_ops[i]

        # Reverse scan: find the first poorly-aligned suffix window.
        index_end = <int> cigar_l
        j = <int> cigar_l - 1
        cov_acc = 0
        mat_acc = 0
        n_acc = 0
        for i in range(<int> cigar_l - 1, index_begin - 1, -1):
            while j > -1 and cov_acc < ws:
                cov_acc += cov_ops[j]
                mat_acc += mat_ops[j]
                n_acc += n_ops[j]
                j -= 1
            if cov_acc >= ws:
                mat_acc_capped = mat_acc if mat_acc < ws else ws
                if mat_acc_capped == 0:
                    w_rate = 0.0
                else:
                    w_rate = n_acc / float(mat_acc_capped + n_acc)
                stop_idx = j
            else:
                if mat_acc == 0:
                    w_rate = 0.0
                else:
                    w_rate = n_acc / float(mat_acc + n_acc)
                stop_idx = -1
            if w_rate < threshold:
                break
            index_end = stop_idx
            cov_acc -= cov_ops[i]
            mat_acc -= mat_ops[i]
            n_acc -= n_ops[i]

        if index_begin <= cigar_index <= index_end:
            spanning.append(item)
    if not spanning_alignments:
        return spanning, 0
    return spanning, 1 - (len(spanning) / len(spanning_alignments))


cpdef gap_size_upper_bound(AlignedSegment alignment, int cigarindex, int pos_input, int end_input, int length_extend=15, float divergence=0.02, float len_t=3):
    # expand indel using cigar, merges nearby gaps into the SV event
    cdef int pos, end, l, opp, extent_left, extent_right, candidate_type, candidate_len, len_input, i, dist, dist_thresh
    cdef uint32_t cigar_value
    cdef uint32_t cigar_l = alignment._delegate.core.n_cigar
    cdef uint32_t *cigar_p = bam_get_cigar(alignment._delegate)

    pos = pos_input
    end = end_input
    extent_left = pos
    extent_right = end
    cigar_value = cigar_p[cigarindex]
    candidate_type = cigar_value & BAM_CIGAR_MASK
    candidate_len = cigar_value >> BAM_CIGAR_SHIFT
    len_input = candidate_len
    dist_thresh = min(<int> (len_input * len_t), <int> (<float> log2_32(1 + candidate_len) / divergence))
    i = cigarindex + 1
    dist = 0

    while i < cigar_l:  # extent rightwards
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & BAM_CIGAR_MASK
        l = <int> cigar_value >> BAM_CIGAR_SHIFT
        if opp == CSOFT_CLIP or opp == CHARD_CLIP:
            break
        if opp == CMATCH or opp == CEQUAL or opp == CDIFF:
            dist += l
            end += l
            if dist > dist_thresh:
                break
        elif l > length_extend and (opp == CDEL or opp == CINS):
            if candidate_type == CDEL:
                if opp == CDEL:
                    candidate_len += l
                    end += l
                    extent_right = end
                else:
                    if l >= candidate_len: break
                    candidate_len -= l
            else:
                if opp == CDEL:
                    if l >= candidate_len: break
                    candidate_len -= l
                    end -= l
                else:
                    candidate_len += l
                    extent_right = end
            dist_thresh = min(<int> (l * len_t), <int> (<float> log2_32(1 + candidate_len) / divergence))
            dist = 0
        i += 1

    i = cigarindex - 1
    dist = 0
    while i > -1:  # extend leftwards
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & BAM_CIGAR_MASK
        l = <int> cigar_value >> BAM_CIGAR_SHIFT
        if opp == CSOFT_CLIP or opp == CHARD_CLIP:
            break
        if opp == CMATCH or opp == CEQUAL or opp == CDIFF:
            dist += l
            pos -= l
            if dist > dist_thresh:
                break
        elif l > length_extend and (opp == CDEL or opp == CINS):
            if candidate_type == CDEL:
                if opp == CDEL:
                    candidate_len += l
                    pos -= l
                    extent_left = pos
                else:
                    if l >= candidate_len: break
                    candidate_len -= l
            else:
                if opp == CDEL:
                    if l >= candidate_len: break
                    candidate_len -= l
                    pos -= l
                else:
                    candidate_len += l
                    extent_left = pos
            dist_thresh = min(<int> (l * len_t), <int>( <float>log2_32(1 + candidate_len) / divergence) )
            dist = 0
        i -= 1

    if extent_right > extent_left and (extent_left != pos_input or extent_right != end_input):  # has been extended
        # middle = extent_left + <int>((extent_right - extent_left) / 2)
        # if candidate_type == CDEL:
        #     pos = middle - int(candidate_len * 0.5)
        #     end = middle + int(candidate_len * 0.5)
        # else:
        #     pos = pos_input
        #     end = end_input
        # return candidate_len, pos, end
        return candidate_len, extent_left, extent_right

    return len_input, pos_input, end_input


def sample_level_density(potential, regions, max_dist=50):
    tmp_list = defaultdict(list)
    cdef EventResult_t ei
    cdef int start, stop, idx
    cdef bint inA, inB
    cdef list region_flags = []

    for idx in range(len(potential)):
        ei = potential[idx]
        # Only find density for non-region calls, otherwise too dense to be meaningful
        inA = intersecter(regions, ei.chrA, ei.posA, ei.posA + 1)
        inB = intersecter(regions, ei.chrB, ei.posB, ei.posB + 1)
        if not inA:
            tmp_list[ei.chrA].append((ei.posA - max_dist, ei.posA + max_dist, idx))
        if not inB:
            tmp_list[ei.chrB].append((ei.posB - max_dist, ei.posB + max_dist, idx))
        region_flags.append((inA, inB))

    si_sets = {}

    for k, v in tmp_list.items():
        iset = IntervalMap()
        for start, stop, idx in v:
            iset.add(start, stop, idx)
        iset.build()
        si_sets[k] = iset

    cdef int vv
    for idx in range(len(potential)):
        ei = potential[idx]
        neighbors = 0.
        count = 0.
        if ei.chrA == ei.chrB and abs(ei.posB - ei.posA) < 2:
            expected = 2
        else:
            expected = 1
        inA, inB = region_flags[idx]
        if not inA:
            vv = si_sets[ei.chrA].count(ei.posA, ei.posA + 1)
            neighbors += vv - expected
            count += 1
        if not inB:
            vv = si_sets[ei.chrB].count(ei.posB, ei.posB + 1)
            neighbors += vv - expected
            count += 1
        neighbors_10kb = 0.
        count_10kb = 0
        large_itv = merge_intervals(((ei.chrA, ei.posA, ei.posA + 1), (ei.chrB, ei.posB, ei.posB + 1)), pad=10000)
        for c, s, e in large_itv:
            if not intersecter(regions, c, s, e):
                vv = si_sets[c].count(s, e)
                neighbors_10kb += vv - len(large_itv)
                count_10kb += 1
        if neighbors < 0:
            neighbors = 0
        if count > 0:
            ei.neigh = neighbors / count
        else:
            ei.neigh = 0
        if count_10kb > 0:
            ei.neigh10kb = neighbors_10kb / count_10kb
        else:
            ei.neigh10kb = 0
    return potential


cdef bint same_k(int start1, int start2, int n, const unsigned char[:] seq):
    cdef int j
    for j in range(n):
        if seq[start1 + j] != seq[start2 + j]:
            return False
    return True


cpdef dict search_ssr_kc(ori):
    cdef bytes seq_b = ori.upper().encode("ascii")
    cdef const unsigned char[:] string_view = seq_b  # use for indexing
    cdef int str_len = len(seq_b)
    cdef int rep_len = min(7, str_len)
    cdef int i = 0
    cdef int t, start, count, mm, good_i, successive_bad, finish
    cdef int n_ref_repeat_bases = 0
    cdef int n_expansion = 0
    cdef int stride = 0
    expansion_seq = ""
    cdef unsigned char starting_base
    cdef int starting_kmer_idx
    cdef int ns, d, shift_ev, shift_total, prev_start
    cdef bint shifted, prev_exact
    cdef list units
    cdef int j, k, span_end, rs, re_, best_size
    cdef bint cur_lower
    cdef list runs
    while i < str_len:
        if seq_b[i] == 78:  # ord('N')
            i += 1
            continue
        for t in range(1, rep_len):
            start = i
            if start + t >= str_len:
                break
            starting_base = string_view[i]
            starting_kmer_idx = start
            count = 1
            mm = 0
            good_i = 0
            successive_bad = 0
            finish = 0
            shift_ev = 0
            shift_total = 0
            prev_exact = True
            units = [i]
            while start + t < str_len and (starting_base == string_view[start+t] or mm < 2):
                start += t
                if start + t > str_len:
                    break
                if not same_k(start, starting_kmer_idx, t, string_view):
                    # Tolerate small indels in the repeat: if the exact-strided
                    # kmer does not match, look for the unit at a small phase
                    # shift. Bounded: not in the first two units, shifts cannot
                    # chain without an intervening exact match, at most 2 shift
                    # events, 3bp total drift, and never across an N
                    shifted = False
                    if shift_ev < 2 and prev_exact and count >= 2:
                        prev_start = start - t
                        for d in (-2, -1, 1, 2):
                            ns = start + d
                            if ns <= prev_start or ns + t > str_len:
                                continue
                            if shift_total + d > 3 or shift_total + d < -3:
                                continue
                            g0 = prev_start + 1
                            g1 = ns + t
                            if g0 > g1:
                                g0, g1 = ns + 1, prev_start + t
                            if seq_b.find(b"N", g0, g1) != -1:
                                continue
                            if same_k(ns, starting_kmer_idx, t, string_view):
                                shifted = True
                                shift_ev += 1
                                shift_total += d
                                start = ns
                                break
                    if shifted:
                        good_i = start
                        finish = good_i
                        successive_bad = 0
                        count += 1
                        units.append(start)
                        prev_exact = False
                        continue
                    successive_bad += 1
                    mm += 1
                    if successive_bad > 1 or count < 2:
                        finish = good_i
                        break
                else:
                    good_i = start
                    finish = good_i
                    successive_bad = 0
                    count += 1
                    units.append(start)
                    prev_exact = True
            if count >= 3 and (finish - i) + t > 10:
                span_end = finish + t
                # Classify the repeat span at base resolution using the
                # lowercase/uppercase convention of the consensus (lowercase =
                # clipped/inserted sequence). A lowercase tract adjacent to an
                # uppercase tract is a repeat expansion at the break site
                runs = []
                rs = i
                cur_lower = ori[i].islower()
                for j in range(i + 1, span_end):
                    if ori[j].islower() != cur_lower:
                        runs.append((rs, j, cur_lower))
                        rs = j
                        cur_lower = not cur_lower
                runs.append((rs, span_end, cur_lower))
                best_size = 0
                for k in range(len(runs)):
                    rs, re_, cur_lower = runs[k]
                    if cur_lower and re_ - rs >= 10:
                        if (k > 0 and not runs[k - 1][2]) or (k < len(runs) - 1 and not runs[k + 1][2]):
                            if re_ - rs > best_size:
                                best_size = re_ - rs
                                n_expansion = re_ - rs
                                expansion_seq = ori[rs:re_]
                                stride = t
                for rs, re_, cur_lower in runs:
                    if not cur_lower:
                        n_ref_repeat_bases += re_ - rs
                i = span_end
        i += 1
    return {"n_expansion": n_expansion, "stride": stride, "exp_seq": expansion_seq, "ref_poly_bases": n_ref_repeat_bases}


def find_repeat_expansions(events, insert_stdev):
    cdef EventResult_t e
    for e in events:
        e.n_expansion = 0
        e.stride = 0
        e.exp_seq = ""
        e.ref_poly_bases = 0
        if e.contig:
            r = search_ssr_kc(e.contig)
            e.n_expansion = r["n_expansion"]
            e.stride = r["stride"]
            e.exp_seq = r["exp_seq"]
            e.ref_poly_bases += r["ref_poly_bases"]
        if e.contig2:
            r = search_ssr_kc(e.contig2)
            e.ref_poly_bases += r["ref_poly_bases"]
            if e.n_expansion < r["n_expansion"]:
                e.n_expansion = r["n_expansion"]
                e.stride = r["stride"]
                e.exp_seq = r["exp_seq"]
    return events
