#cython: language_level=3, boundscheck=False
from __future__ import absolute_import
from collections import deque, defaultdict
from typing import Dict, DefaultDict, List, Set, Tuple
import logging

import superintervals
from dysgu import io_funcs
import numpy as np
cimport numpy as np
from scipy.optimize import nnls
import pysam
DTYPE = np.float64
ctypedef np.float_t DTYPE_t
from dysgu.map_set_utils cimport CoverageTrack
from dysgu.map_set_utils cimport TranscriptData
from dysgu.map_set_utils import merge_intervals, echo
from dysgu.io_funcs import intersecter

from libc.stdint cimport uint32_t
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_cigar
from libcpp.string cimport string
from libcpp.vector cimport vector


def index_stats(f, rl=None):
    if rl is None:
        rl = []
        for idx, a in enumerate(f.fetch()):
            if not a.flag & 3840 and a.cigartuples is not None:
                rl.append(a.infer_read_length())
            if idx > 10000:
                break
        rl = np.median(rl)
    s = f.get_index_statistics()
    total = 0
    ref_len = 0
    for i in s:
        total += i.mapped
        ref_len += f.get_reference_length(i.contig)
    cov = (total*rl/ref_len)*0.98
    return cov, rl


def auto_max_cov(mc, bname):
    if mc == "-1":
        mc = 1e6
    elif mc == "auto":
        if bname == "-":
            raise NotImplementedError("Not possible to use max-cov == 'auto' with stdin")
        aln_f = pysam.AlignmentFile(bname)
        if aln_f.is_cram:
            raise NotImplementedError("Not possible to use index_stats on a cram file")
        mc, rl = index_stats(aln_f)
        mc = round(mc) * 6
        if mc < 5:
            logging.critical("Max-cov estimated as < 5? Set manually to proceed")
        logging.info("Auto max-cov estimated {}x".format(mc))
    else:
        try:
            mc = int(mc)
        except:
            raise ValueError("Max-cov value not understood {}".format(mc))
    return mc


def median(L):
    if len(L) % 2 == 1:
        return L[int(len(L)/2)]
    mid = int(len(L) / 2) - 1
    return (L[mid] + L[mid+1]) / 2.0


def unscaled_upper_mad(xs):
    """Return a tuple consisting of the median of xs followed by the
    unscaled median absolute deviation of the values in xs that lie
    above the median.
    """
    xs.sort()
    med = median(xs)
    umad = median([x - med for x in xs if x > med])
    return med, umad


def mean_std(L):
    s = sum(L)
    mean = np.median(L)
    sq_sum = 0.0
    for v in L:
        sq_sum += (v - mean)**2.0
    var = sq_sum / float(len(L))
    return mean, var**0.5


def get_insert_params(L, mads=8):  # default for lumpy is 10
    c = len(L)
    med, umad = unscaled_upper_mad(L)
    upper_cutoff = med + mads * umad
    L = [v for v in L if v < upper_cutoff]
    new_len = len(L)
    removed = c - new_len
    logging.info("Calculating insert size. Removed {} outliers with insert size >= {}".format(removed, upper_cutoff))
    mean, stdev = mean_std(L)
    return mean, stdev


cdef class GenomeScanner:
    def __init__(self, inputbam, int mapq_threshold, int max_cov, include_regions, read_threads, buffer_size, regions_only, stdin,
                 clip_length=30, min_within_size=30, cov_track_path=None, paired_end=True, bam_iter=None):
        self.input_bam = inputbam
        self.mapq_threshold = mapq_threshold
        self.max_cov = max_cov
        self.include_regions = include_regions
        self.overlap_regions = io_funcs.overlap_regions(include_regions, int_chroms=True, infile=inputbam)
        self.regions_only = regions_only
        self.clip_length = clip_length
        self.min_within_size = min_within_size
        self.staged_reads = deque([])
        self.procs = read_threads
        self.buff_size = buffer_size
        self.current_bin = []
        self.current_cov_array = None
        self.reads_dropped = 0
        self.depth_d = {}
        self.cov_track_path = cov_track_path
        self.first = 1  # Not possible to get first position in a target fetch region, so buffer first read instead
        self.read_buffer = dict()
        self.approx_read_length = -1
        self.last_tell = 0
        self.no_tell = True if stdin else False
        self.paired_end = paired_end
        self.bam_iter = bam_iter
        self.cpp_cov_track = CoverageTrack()
        self.current_tid = -1

    def _get_reads(self):
        # Two options, reads are collected from whole genome, or from target regions only
        cdef int index_start, opp, length, chrom_length, pos, i
        cdef bytes out_path
        cdef bint good_read
        cdef int mq_thresh = self.mapq_threshold
        cdef AlignedSegment aln
        cdef uint32_t cigar_value
        cdef uint32_t cigar_l
        cdef uint32_t *cigar_p

        if not self.include_regions or not self.regions_only:
            # Some reads may have been staged from getting read_length (if file is streamed from stdin)
            while len(self.staged_reads) > 0:
                yield self.staged_reads.popleft()
            tell = 0 if self.no_tell else self.input_bam.tell()
            # when 'run' command is invoked, run this block. cov track already exists from find-reads
            if self.cov_track_path is None and self.bam_iter is None:
                for aln in self.input_bam:
                    cigar_l = aln._delegate.core.n_cigar

                    if aln.flag & 1284 or aln.mapq < mq_thresh or cigar_l == 0:  # not primary, duplicate or unmapped?
                        continue
                    yield aln, tell
                    tell = 0 if self.no_tell else self.input_bam.tell()
                    # yield aln, tell

            # when 'call' command was run only. generate coverage track here
            else:
                self.cpp_cov_track.set_max_cov(self.max_cov)
                if self.bam_iter is not None:
                    f_iter = self.bam_iter
                else:
                    f_iter = self.input_bam.fetch(until_eof=True)
                for aln in f_iter:

                    cigar_l = aln._delegate.core.n_cigar
                    if cigar_l == 0 or aln.flag & 1284:
                        # tell = 0 if self.no_tell else self.input_bam.tell()
                        continue

                    if aln.rname != self.current_tid:
                        if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                            out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                            self.cpp_cov_track.write_track(out_path)
                        chrom_length = self.input_bam.get_reference_length(self.input_bam.get_reference_name(aln.rname))
                        self.current_tid = aln.rname
                        self.cpp_cov_track.set_cov_array(chrom_length)
                    # Add to coverage track here
                    index_start = 0
                    pos = aln.pos
                    good_read = False
                    if aln.flag & 2048:
                        good_read = True
                    elif not aln.flag & 2:
                        good_read = True
                    cigar_p = bam_get_cigar(aln._delegate)
                    for i in range(cigar_l):
                        cigar_value = cigar_p[i]
                        opp = <int> cigar_value & 15
                        length = <int> cigar_value >> 4
                        if opp == 4 and length >= self.clip_length:
                            good_read = True
                        elif opp == 2 or opp == 3:
                            index_start += length
                            if length >= self.min_within_size:
                                good_read = True
                        elif opp == 1 and length >= self.min_within_size:
                            good_read = True
                        elif opp == 0 or opp == 7 or opp == 8:
                            self.cpp_cov_track.add(pos + index_start, pos + index_start + length)
                            index_start += length

                    if aln.mapq < mq_thresh or not good_read or not self.cpp_cov_track.cov_val_good(self.current_tid, aln.rname, pos):
                        continue
                    yield aln, tell
                    tell = 0 if self.no_tell else self.input_bam.tell()
                    # yield aln, tell

                if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                    out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                    self.cpp_cov_track.write_track(out_path)

        else:
            logging.info("Searching --regions and collecting reads from mate-pair locations")
            # This block is for when --regions-only is set on the call command
            self.cpp_cov_track.set_max_cov(self.max_cov)
            # Reads must be fed into graph in sorted order, find regions of interest first
            intervals_to_check = []  # Containing include_regions, and also mate pairs
            pad = 1000
            regions = [j.strip().split("\t")[:3] for j in open(self.include_regions, "r") if j[0] != "#"]
            for c, s, e in regions:
                intervals_to_check.append((c, int(s), int(e)))
            for c, s, e in regions:
                for a in self.input_bam.fetch(c, int(s), int(e)):
                    if not a.flag & 1800:
                        p1 = a.pnext - pad
                        c2 = self.input_bam.get_reference_name(a.rnext)
                        if c2:
                            intervals_to_check.append((c2, 0 if p1 < 0 else p1, a.pnext + pad))
                        if a.has_tag("SA"):
                            sa = a.get_tag("SA").split(";")
                            for v in sa:
                                try:
                                    chrom2, pos2, _ = v.split(",", 2)
                                except ValueError:  # sometimes trailing empty string in list ""
                                    break
                                pos2 = int(pos2) - pad
                                intervals_to_check.append((chrom2, 0 if pos2 < 0 else pos2, pos2 + pad))

            itv = merge_intervals(intervals_to_check)
            seen_reads = set([])  # Avoid reading same alignment twice, shouldn't happen anyway
            for c, s, e in itv:
                tell = -1  # buffer first read because tell will likely be wrong
                for aln in self.input_bam.fetch(c, int(s), int(e)):
                    name = aln.qname.__hash__(), aln.flag, aln.pos
                    if name in seen_reads:
                        continue

                    cigar_l = aln._delegate.core.n_cigar
                    if aln.flag & 1284 or aln.mapq < mq_thresh or cigar_l == 0:
                        continue
                    if aln.rname != self.current_tid:
                        if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                            out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                            self.cpp_cov_track.write_track(out_path)
                        chrom_length = self.input_bam.get_reference_length(self.input_bam.get_reference_name(aln.rname))
                        self.current_tid = aln.rname
                        self.cpp_cov_track.set_cov_array(chrom_length)
                    index_start = 0
                    pos = aln.pos
                    good_read = False
                    if aln.flag & 2048:
                        good_read = True
                    elif not aln.flag & 2:
                        good_read = True
                    cigar_l = aln._delegate.core.n_cigar
                    cigar_p = bam_get_cigar(aln._delegate)
                    for i in range(cigar_l):
                        cigar_value = cigar_p[i]
                        opp = <int> cigar_value & 15
                        length = <int> cigar_value >> 4
                        if opp == 4 and length >= self.clip_length:
                            good_read = True
                        elif opp == 2 or opp == 3:
                            index_start += length
                            if length >= self.min_within_size:
                                good_read = True
                        elif opp == 1 and length >= self.min_within_size:
                            good_read = True
                        elif opp == 0 or opp == 7 or opp == 8:
                            self.cpp_cov_track.add(pos + index_start, pos + index_start + length)
                            index_start += length
                    seen_reads.add(name)
                    if not good_read:
                        continue
                    if not self.cpp_cov_track.cov_val_good(self.current_tid, aln.rname, pos):
                        continue
                    yield aln, tell
                    tell = 0 if self.no_tell else self.input_bam.tell()
                    # yield aln, tell

                if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                    out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                    self.cpp_cov_track.write_track(out_path)


    def get_read_properties(self, int max_tlen, int insert_median, int insert_stdev, int read_len, ibam=None, find_divergence=False):
        # This is invoked first to scan the first part of the file for the insert size metrics,
        # or open and process the --ibam alignment file
        if not find_divergence and read_len != -1:
            logging.info(f"Read length {read_len}, insert_median {insert_median}, insert stdev {insert_stdev}")
            self.approx_read_length = read_len
            return insert_median, insert_stdev
        approx_read_length_l = []
        inserts = []
        divergence = []
        cdef int required = 97
        restricted = 3484
        cdef int flag_mask = required | restricted
        cdef int c = 0
        tell = 0
        cdef int flag, tlen
        cdef float approx_read_length
        if ibam is not None:
            file_iter = ibam
        else:
            file_iter = self.input_bam
        cdef bint any_paired_end = False
        cdef AlignedSegment a
        cdef uint32_t cigar_value
        cdef uint32_t cigar_l
        cdef uint32_t *cigar_p
        cdef int i
        cdef int non_match_count, matched_bases
        
        # first check SO tag in HD
        if "HD" in file_iter.header:
            hd = file_iter.header["HD"]
            if "SO" in hd:
                if hd["SO"] == "unsorted":
                    raise ValueError("Input bam file must be sorted")

        for a in file_iter:
            if ibam is None:
                cigar_l = a._delegate.core.n_cigar
                if a.flag & 1284 or a.mapq < self.mapq_threshold or cigar_l == 0:
                    continue
                tell = 0 if self.no_tell else self.input_bam.tell()
                if self.no_tell:
                    self.staged_reads.append((a, tell))
                    # self._add_to_bin_buffer(a, tell)
                    if a.rname != self.current_tid:
                        if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                            out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path.outpath, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                            self.cpp_cov_track.write_track(out_path)
                        chrom_length = self.input_bam.get_reference_length(self.input_bam.get_reference_name(a.rname))
                        self.current_tid = a.rname
                        self.cpp_cov_track.set_cov_array(chrom_length)
                    index_start = 0
                    pos = a.pos
                    cigar_l = a._delegate.core.n_cigar
                    cigar_p = bam_get_cigar(a._delegate)
                    for i in range(cigar_l):
                        cigar_value = cigar_p[i]
                        opp = <int> cigar_value & 15
                        length = <int> cigar_value >> 4
                        if opp == 2 or opp == 3:
                            index_start += length
                        elif opp == 0 or opp == 7 or opp == 8:
                            end = index_start + length
                            self.cpp_cov_track.add(pos + index_start, pos + end)
                            index_start += length
            if find_divergence:
                if len(divergence) < 20000:
                    non_match_count = 0
                    matched_bases = 0
                    cigar_l = a._delegate.core.n_cigar
                    cigar_p = bam_get_cigar(a._delegate)
                    for i in range(cigar_l):
                        cigar_value = cigar_p[i]
                        opp = <int> cigar_value & 15
                        length = <int> cigar_value >> 4
                        if opp == 0 or opp == 7:
                            matched_bases += length
                        elif opp == 1 or opp == 2 or opp == 8:
                            non_match_count += 1
                    divergence.append(non_match_count / (non_match_count + matched_bases))
                else:
                    break
            else:
                if len(approx_read_length_l) < 200000:
                    flag = a.flag
                    if a.seq is not None:
                        rl = a.infer_read_length()
                        if rl:
                            approx_read_length_l.append(rl)
                            if a.rname == a.rnext and flag & flag_mask == required and a.tlen >= 0:
                                inserts.append(a.tlen)
                else:
                    break
                if c > 20000000:
                    logging.critical("Cant infer read properties after 10 million reads, please set manually")
                    return -1, -1
            c += 1

        if find_divergence and len(divergence) == 0:
            logging.critical("Cant infer read divergence, no reads?")
            return -1, -1
        elif not find_divergence and len(approx_read_length_l) == 0:
            logging.critical("Cant infer read insert and size, no reads?")
            return -1, -1
        if find_divergence:
            div = divergence[-max(0, len(divergence) - 500):]
            mean = np.mean(divergence)
            divergence_upper_bound = min(1, max(2.5 * mean, mean + (2.5 * np.std(divergence))))
            logging.info(f"Inferred sequence divergence upper bound {round(divergence_upper_bound,4)}")
            if ibam is None:
                self.last_tell = tell
                if not self.no_tell:
                    self.input_bam.reset()
            return divergence_upper_bound, 0
        else:
            approx_read_length = int(np.median(approx_read_length_l))
            self.approx_read_length = approx_read_length
            assume_single = False
            if len(inserts) <= 100 and insert_median == -1:
                insert_median = 300
                insert_stdev = 150
                assume_single = True
            if insert_median == -1:
                insert_median, insert_stdev = get_insert_params(inserts)
            if not assume_single:
                logging.info(f"Inferred read length: {approx_read_length}, insert median: {insert_median}, insert stdev: {insert_stdev}")
            else:
                logging.info(f"Inferred read length: {approx_read_length}")
            if ibam is None:
                self.last_tell = tell
                if not self.no_tell:
                    self.input_bam.reset()
            return insert_median, insert_stdev

    def iter_genome(self):
        # Read the rest of the genome, reads are sent in blocks
        cdef int total_reads = 0
        for staged in self._get_reads():
            total_reads += 1  # len(staged)
            yield staged
        if total_reads == 0:
            logging.critical("No reads found, finishing")
            return
        logging.info(f"Total input reads: {total_reads}")

    def add_to_buffer(self, r, n1, tell):
        if self.first == 1 or tell == -1:
            # Not possible to get tell position from first read in region, so put into buffer instead
            self.read_buffer[n1] = r
            self.first = 0
        elif len(self.read_buffer) < self.buff_size:
            self.read_buffer[n1] = r
        elif self.no_tell:
            raise BufferError("Read buffer has overflowed, increase --buffer-size")


cdef float add_coverage(int start, int end, DTYPE_t[:] chrom_depth) nogil:
    cdef float fs = start / 100
    cdef float fe = end / 100
    cdef int bin_start = <int> fs
    cdef int bin_end = <int> fe
    if bin_start < 0:
        bin_start = 0
    if bin_end > <int> len(chrom_depth) - 1:
        bin_end = <int> len(chrom_depth) - 1
    cdef float ol_start = <float> bin_start + 1 - fs
    chrom_depth[bin_start] += ol_start
    cdef float ol_end = 0
    if bin_start != bin_end:
        ol_end = fe - (<float> bin_end)
        chrom_depth[bin_end] += ol_end
        if bin_end - bin_start > 1:
            for i in range(bin_start + 1, bin_end):
                chrom_depth[i] += 1
    return chrom_depth[bin_start]


cpdef calculate_coverage(int start, int end, np.int16_t[:] chrom_depth, int bin_size=10, bint ignore_zero=0):
    cdef float fs = start / bin_size
    cdef float fe = end / bin_size
    start = <int> fs
    end = <int> fe
    if start < 0:
        start = 0

    cdef int len_chrom = <int> len(chrom_depth)
    if end > len_chrom:
        end = len_chrom
    if start >= len_chrom:
        return 0, 0
    if end < 0:
        return 0, 0
    if end < start:
        end = start

    cdef int i
    cdef float total = 0
    cdef float max_cov = 0
    cdef float cov_val = 0
    cdef int n_bins = 0  # number of bins counted in mean
    with nogil:
        if start == end:
            cov_val = chrom_depth[start]
            if (not ignore_zero) or cov_val != 0:
                total = cov_val
                n_bins = 1
                max_cov = cov_val
            else:
                total = 0
                n_bins = 0
                max_cov = 0
        else:
            for i in range(start, end):
                cov_val = chrom_depth[i]
                if cov_val > max_cov:
                    max_cov = cov_val
                if ignore_zero and cov_val == 0:
                    continue
                total += cov_val
                n_bins += 1

    if total == 0 or n_bins == 0:
        return 0, 0

    if start == end:
        return total, max_cov

    return total / n_bins, max_cov



def switch_AB(r):
    chrA, posA, cipos95A, contig2 = r.chrA, r.posA, r.cipos95A, r.contig2
    r.chrA = r.chrB
    r.posA = r.posB
    r.cipos95A = r.cipos95B
    r.chrB = chrA
    r.posB = posA
    r.cipos95B = cipos95A
    r.contig2 = r.contig
    r.contig = contig2


def get_raw_coverage_information(events, regions, regions_depth):
    new_events = []
    for r in events:
        ar = False
        if intersecter(regions, r.chrA, r.posA, r.posA + 1):
            ar = True
        br = False
        if intersecter(regions, r.chrB, r.posB, r.posB + 1):
            br = True
        kind = "extra-regional"

        if not ar and not br:
            if r.chrA == r.chrB and r.posA > r.posB:  # Put non-region first
                switch_AB(r)
        if (br and not ar) or (not br and ar):
            kind = "hemi-regional"
            if not br and ar:
                switch_AB(r)

        if ar and br:
            if r.chrA == r.chrB:
                rA = regions[r.chrA].search_values(r.posA, r.posA + 1)
                rB = regions[r.chrB].search_values(r.posB, r.posB + 1)
                if rA and rB and rA[0] == rB[0]:
                    kind = "intra_regional"
                    if r.posA > r.posB:
                        switch = True
                else:
                    kind = "inter-regional"
                    if r.chrA != sorted([r.chrA, r.chrB])[0]:
                        switch = True
            else:
                kind = "inter-regional"

        max_depth = 0
        if kind == "hemi-regional":
            chrom_i = r.chrA
            if chrom_i in regions_depth.chrom_cov_arrays:
                reads_10kb, max_depth = calculate_coverage(r.posA - 10000, r.posA + 10000, regions_depth.chrom_cov_arrays[chrom_i], ignore_zero=0)
                reads_10kb = round(reads_10kb, 3)
            else:
                reads_10kb = 0
        else:
            chrom_i = r.chrA
            if chrom_i in regions_depth.chrom_cov_arrays:
                reads_10kb_left, max_depth = calculate_coverage(r.posA - 10000, r.posA + 10000, regions_depth.chrom_cov_arrays[chrom_i], ignore_zero=0)
                reads_10kb_left = round(reads_10kb_left, 3)
            else:
                reads_10kb_left = 0
            chrom_i = r.chrB
            if chrom_i in regions_depth.chrom_cov_arrays:
                reads_10kb_right, max_depth = calculate_coverage(r.posB - 10000, r.posB + 10000, regions_depth.chrom_cov_arrays[chrom_i], ignore_zero=0)
                reads_10kb_right = round(reads_10kb_right, 3)
            else:
                reads_10kb_right = 0
            if reads_10kb_left > reads_10kb_right:
                reads_10kb = reads_10kb_left
            else:
                reads_10kb = reads_10kb_right

        r.kind = kind
        r.raw_reads_10kb = reads_10kb
        if r.chrA != r.chrB:
            r.svlen = 1000000
        r.mcov = max_depth

        r.a_freq = r.a_freq / (reads_10kb + 1e-5)
        r.a_freq = round(max(0, min(r.a_freq, 1.0)), 3)
        new_events.append(r)
    return new_events


def transcriptome_effective_depth(
    chrom: str,
    pos: int,
    chrom_depth,
    tr_intervals,
    blocks_by_tx,
    bin_size: int = 10,
    flank: int = 10_000,        # only use exon bins near breakpoint
    min_bins: int = 30,       # fallback if too few bins
):
    # Find transcripts overlapping breakpoint
    txs = tr_intervals.get(chrom, None)
    if txs is None:
        return 0.0, 0.0

    hits = txs.search_values(pos, pos + 1)
    if not hits:
        return 0.0, 0.0

    # Collect bin indices from exon blocks near breakpoint
    lo = pos - flank
    hi = pos + flank
    n = int(len(chrom_depth))
    bins = set()
    for tx in set(hits):
        for s, e in blocks_by_tx.get(tx, ()):
            if e < lo or s > hi:
                continue
            bs = s // bin_size
            be = (e + bin_size - 1) // bin_size  # ceil
            if bs < 0:
                bs = 0
            if be > n:
                be = n
            for b in range(bs, be):
                bins.add(b)

    if len(bins) < min_bins:
        return 0.0, 0.0

    depths = np.fromiter((chrom_depth[b] for b in bins), dtype=np.int32)
    nz = depths[depths > 0]
    if nz.size == 0:
        return 0.0, 0.0

    # Use bins above noise floor
    thr = float(np.percentile(nz, 1))
    if thr < 1.0:
        thr = 1.0
    kept = nz[nz >= thr]
    if kept.size == 0:
        return float(nz.mean()), float(nz.max())

    return float(kept.mean()), float(kept.max())


def load_transcript_blocks(transcript_blocks_path: str):
    tr_intervals = defaultdict(lambda: superintervals.IntervalMap())
    blocks_by_tx = defaultdict(list)
    tx_strand: Dict[str, str] = {}
    tx_chrom: Dict[str, str] = {}

    with open(transcript_blocks_path, "r") as tb:
        for line in tb:
            if not line or line[0] == "#":
                continue
            l = line.rstrip("\n").split("\t")
            if len(l) < 4:
                continue
            chrom = l[0]
            start = int(l[1])
            end = int(l[2])
            tx = l[3]
            strand = l[5] if len(l) > 5 else "+"

            blocks_by_tx[tx].append((start, end))
            tx_strand[tx] = strand
            tx_chrom[tx] = chrom

    # Add transcript span intervals per chromosome
    for tx, blocks in blocks_by_tx.items():
        if not blocks:
            continue
        chrom = tx_chrom[tx]
        min_start = min(s for s, _ in blocks)
        max_end = max(e for _, e in blocks)
        tr_intervals[chrom].add(min_start, max_end, tx)

    for imap in tr_intervals.values():
        imap.build()

    return tr_intervals, blocks_by_tx, tx_strand, tx_chrom


def fallback_genomic(chrom, pos, regions_depth):
    if chrom in regions_depth.chrom_cov_arrays:
        mean_cov, mx = calculate_coverage(pos - 10000, pos + 10000, regions_depth.chrom_cov_arrays[chrom], ignore_zero=1)
        return float(mean_cov), float(mx)
    return 0.0, 0.0


def depth_at(chrom, pos, regions_depth, tr_intervals, blocks_by_tx):
    if chrom not in regions_depth.chrom_cov_arrays:
        return 0.0, 0.0
    mean_cov, mx = transcriptome_effective_depth(
        chrom, pos,
        regions_depth.chrom_cov_arrays[chrom],
        tr_intervals, blocks_by_tx,
        bin_size=10,
        flank=10_000,
        min_bins=10
    )
    if mean_cov == 0.0:
        return fallback_genomic(chrom, pos, regions_depth)
    return mean_cov, mx


def common_transcript_overlap(tr_intervals, chrA, posA, chrB, posB) -> bool:
    if chrA != chrB:
        return False
    imap = tr_intervals.get(chrA)
    if imap is None:
        return False
    hitsA = imap.search_values(posA, posA + 1) or []
    hitsB = imap.search_values(posB, posB + 1) or []
    return bool(set(hitsA) & set(hitsB))


def in_any_transcript(tr_intervals, chrom, pos) -> bool:
    if chrom not in tr_intervals[chrom]:
        return False
    return tr_intervals[chrom].has_overlaps(pos, pos + 1)


def collect_directional_nonzero_bins(cov_array, pos, direction, max_keep, scan_limit_bins=400) -> list:
    # Collect up to max_keep NONZERO bins, scanning up to scan_limit_bins bins.
    n = int(len(cov_array))
    b0 = int(pos // 10)  # 10=bin_size
    depths = []
    if direction == "left":
        i = b0
        steps = 0
        while i >= 0 and steps < scan_limit_bins and len(depths) < max_keep:
            d = float(cov_array[i])
            if d > 0:
                depths.append(d)
            i -= 1
            steps += 1
    elif direction == "right":
        i = b0
        steps = 0
        while i < n and steps < scan_limit_bins and len(depths) < max_keep:
            d = float(cov_array[i])
            if d > 0:
                depths.append(d)
            i += 1
            steps += 1
    else:
        raise ValueError("direction must be 'left' or 'right'")
    return depths


def metrics_from_depths(depths: list, use_noise_floor=False):
    if not depths:
        return {"thr": None, "n_nonzero": 0, "n_kept": 0, "mean": 0.0, "cv": 0.0}
    arr = np.asarray(depths, dtype=np.float64)
    if use_noise_floor:
        thr = float(np.percentile(arr, 2.))
        if thr > 0.5:
            kept = arr[arr >= thr]
        else:
            kept = arr
    else:
        thr = 1.0
        kept = arr
    n_kept = int(kept.size)
    if n_kept == 0:
        return {"thr": float(thr), "n_nonzero": int(arr.size), "n_kept": 0, "mean": 0.0, "cv": 0.0}

    mean = float(np.mean(kept))
    if kept.size > 1 and mean > 0:
        cv = float(np.std(kept) / mean)
    else:
        cv = 0.0
    return {"thr": float(thr), "n_nonzero": int(arr.size), "n_kept": n_kept, "mean": mean, "cv": cv}


# ------------------
# Entry for RNAseq
# ------------------
def get_raw_coverage_information_transcriptome(events, regions, regions_depth, transcript_blocks_path):
    assert transcript_blocks_path

    max_bins_in_tx: int = 50       # if breakpoint overlaps any transcript, collect up to this many bins
    max_bins_out_tx: int = 15      # otherwise collect up to this many bins

    tr_intervals, blocks_by_tx, tx_strand, tx_chrom = load_transcript_blocks(transcript_blocks_path)

    new_events = []
    for r in events:

        reads_left, max_left = depth_at(r.chrA, r.posA, regions_depth, tr_intervals, blocks_by_tx)
        reads_right, max_right = depth_at(r.chrB, r.posB, regions_depth, tr_intervals, blocks_by_tx)

        reads_10kb = reads_left if reads_left > reads_right else reads_right
        max_depth = max(max_left, max_right)

        r.kind = "extra-regional"
        r.raw_reads_10kb = round(reads_10kb, 3)
        if r.chrA != r.chrB:
            r.svlen = 1000000
        r.mcov = max_depth
            
        if r.chrA not in regions_depth.chrom_cov_arrays or r.chrB not in regions_depth.chrom_cov_arrays:
            new_events.append(r)
            continue

        covA = regions_depth.chrom_cov_arrays[r.chrA]
        covB = regions_depth.chrom_cov_arrays[r.chrB]

        # Choose max bins based on whether each breakpoint lies in any transcript
        max_left = max_bins_in_tx if in_any_transcript(tr_intervals, r.chrA, r.posA) else max_bins_out_tx
        max_right = max_bins_in_tx if in_any_transcript(tr_intervals, r.chrB, r.posB) else max_bins_out_tx

        direction_A = "left" if r.join_type[0] == '3' else "right"
        direction_B = "left" if r.join_type[-1] == '3' else "right"
        left_depths = collect_directional_nonzero_bins(covA, r.posA, direction_A, max_keep=max_left)
        right_depths = collect_directional_nonzero_bins(covB, r.posB, direction_B, max_keep=max_right)

        L = metrics_from_depths(left_depths, use_noise_floor=True)
        R = metrics_from_depths(right_depths, use_noise_floor=True)

        reads_left = L['mean']  # Get directional coverage 
        reads_right = R['mean']

        same_tx = common_transcript_overlap(tr_intervals, r.chrA, r.posA, r.chrB, r.posB)
        if same_tx:
            denom = 0.5 * (reads_left + reads_right)
        else:  # Use the min side
            denom = min(reads_left, reads_right)
        r.a_freq = r.a_freq / (denom + 1e-5)
        r.a_freq = round(max(0, min(r.a_freq, 1.0)), 3)

        # todo this hard filter may need to be a parameter
        ratio = R['mean'] / (L['mean']+1e-6) if (L['mean'] > R['mean']) else L['mean'] / (R['mean']+1e-6)
        keep_event = ratio > 0.01
        if keep_event:
            new_events.append(r)
            # echo(f"KEPT: {r.chrA}:{r.posA}-{r.chrB}:{r.posB} su={r.su} af={r.a_freq}", L, R)
        # else:
            # echo(f"FILTERED: {r.chrA}:{r.posA}-{r.chrB}:{r.posB} su={r.su} af={r.a_freq}", L, R)
        # new_events.append(r)

    return new_events

