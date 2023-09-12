#cython: language_level=3, boundscheck=False
from __future__ import absolute_import
from collections import deque
from dysgu import io_funcs
import numpy as np
cimport numpy as np
import logging
import pysam
DTYPE = np.float64
ctypedef np.float_t DTYPE_t
from dysgu.map_set_utils cimport CoverageTrack
from dysgu.map_set_utils import merge_intervals, echo
from dysgu.io_funcs import intersecter
from libc.stdint cimport uint32_t
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_cigar

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
        self.current_cov = 0
        self.current_chrom = 0
        self.current_pos = 0
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
                    if aln.flag & 1284 or aln.mapq < mq_thresh or aln.cigartuples is None:  # not primary, duplicate or unmapped?
                        continue
                    self._add_to_bin_buffer(aln, tell)
                    tell = 0 if self.no_tell else self.input_bam.tell()
                    while len(self.staged_reads) > 0:
                        yield self.staged_reads.popleft()
            # when 'call' command was run only. generate coverage track here
            else:
                self.cpp_cov_track.set_max_cov(self.max_cov)
                if self.bam_iter is not None:
                    f_iter = self.bam_iter
                else:
                    f_iter = self.input_bam.fetch(until_eof=True)
                for aln in f_iter:
                    # if aln.flag & 1284 or aln.mapq < mq_thresh or aln.cigartuples is None:
                    #     continue
                    cigar_l = aln._delegate.core.n_cigar
                    if cigar_l == 0 or aln.flag & 1284:
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
                    # cigar_l = aln._delegate.core.n_cigar
                    cigar_p = bam_get_cigar(aln._delegate)
                    for i in range(cigar_l):
                        cigar_value = cigar_p[i]
                        opp = <int> cigar_value & 15
                        length = <int> cigar_value >> 4
                        if opp == 4 and length >= self.clip_length:
                            good_read = True
                        elif opp == 2:
                            index_start += length
                            if length >= self.min_within_size:
                                good_read = True
                        elif opp == 1 and length >= self.min_within_size:
                            good_read = True
                        elif opp == 0 or opp == 7 or opp == 8:
                            self.cpp_cov_track.add(pos + index_start, pos + index_start + length)
                            index_start += length
                    if aln.mapq < mq_thresh:
                        continue
                    if not good_read:
                        continue
                    if not self.cpp_cov_track.cov_val_good(self.current_tid, aln.rname, pos):
                        continue
                    self._add_to_bin_buffer(aln, tell)
                    tell = 0 if self.no_tell else self.input_bam.tell()
                    while len(self.staged_reads) > 0:
                        yield self.staged_reads.popleft()
                if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                    out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                    self.cpp_cov_track.write_track(out_path)
            if len(self.current_bin) > 0:
                yield self.current_bin

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
                    if aln.flag & 1284 or aln.mapq < mq_thresh or aln.cigartuples is None:
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
                        elif opp == 2:
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
                    self._add_to_bin_buffer(aln, tell)
                    tell = 0 if self.no_tell else self.input_bam.tell()
                    while len(self.staged_reads) > 0:
                        yield self.staged_reads.popleft()

                if self.current_tid != -1 and self.current_tid <= self.input_bam.nreferences:
                    out_path = "{}/{}.dysgu_chrom.bin".format(self.cov_track_path, self.input_bam.get_reference_name(self.current_tid)).encode("ascii")
                    self.cpp_cov_track.write_track(out_path)
                if len(self.current_bin) > 0:
                    yield self.current_bin


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
        for a in file_iter:
            if ibam is None:
                if a.flag & 1284 or a.mapq < self.mapq_threshold or a.cigartuples is None:
                    continue
                tell = 0 if self.no_tell else self.input_bam.tell()
                if self.no_tell:
                    self._add_to_bin_buffer(a, tell)
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
                        if opp == 2:
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
                logging.info(f"Inferred read length {approx_read_length}, insert median {insert_median}, insert stdev {insert_stdev}")
            else:
                logging.info(f"Inferred read length {approx_read_length}")
            if ibam is None:
                self.last_tell = tell
                if not self.no_tell:
                    self.input_bam.reset()
            return insert_median, insert_stdev

    def iter_genome(self):
        # Read the rest of the genome, reads are sent in blocks
        cdef int total_reads = 0
        for staged in self._get_reads():
            total_reads += len(staged)
            yield staged
        if total_reads == 0:
            logging.critical("No reads found, finishing")
            return
        logging.info(f"Total input reads {total_reads}")

    def add_to_buffer(self, r, n1, tell):
        if self.first == 1 or tell == -1:
            # Not possible to get tell position from first read in region, so put into buffer instead
            self.read_buffer[n1] = r
            self.first = 0
        elif len(self.read_buffer) < self.buff_size:
            self.read_buffer[n1] = r
        elif self.no_tell:
            raise BufferError("Read buffer has overflowed, increase --buffer-size")

    def _add_to_bin_buffer(self, a, tell):
        # Calculates coverage information on fly, drops high coverage regions, buffers reads
        cdef int flag = a.flag
        if flag & 1540 or a.cigartuples is None or a.seq is None:
            return
        cdef int rname = a.rname
        cdef int apos = a.pos
        cdef int bin_pos = int(apos / 100)
        cdef int ref_length
        cdef str reference_name = ""
        cdef int aend = a.reference_end
        cdef float current_coverage
        if self.current_chrom != rname:
            self.current_chrom = rname
        in_roi = False
        if self.overlap_regions:
            in_roi = intersecter(self.overlap_regions, a.rname, apos, apos+1)
        if rname == self.current_chrom and bin_pos == self.current_pos:
            self.current_bin.append((a, tell))
        else:
            if len(self.current_bin) != 0:
                self.staged_reads.append(self.current_bin)
            self.current_chrom = rname
            self.current_pos = bin_pos
            self.current_bin = [(a, tell)]


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


cpdef calculate_coverage(int start, int end, np.int16_t[:] chrom_depth, int bin_size=10):
    cdef float fs = start / bin_size
    cdef float fe = end / bin_size
    start = <int> fs
    end = <int> fe
    if start < 0:
        start = 0
    cdef int len_chrom = <int> len(chrom_depth)
    if end > len_chrom:
        end = len_chrom
    cdef int i
    cdef float total = 0
    cdef float max_cov = 0
    cdef float cov_val = 0
    with nogil:
        if start == end:
            total = chrom_depth[start]
            max_cov = total
        else:
            for i in range(start, end):
                cov_val = chrom_depth[i]
                total += cov_val
                if cov_val > max_cov:
                    max_cov = cov_val
    if total == 0:
        return 0, 0
    if start == end:
        return total, max_cov
    return total / (end - start), max_cov


def get_raw_coverage_information(events, regions, regions_depth, infile, max_cov):
    new_events = []
    for r in events:
        ar = False
        if intersecter(regions, r.chrA, r.posA, r.posA + 1):
            ar = True
        br = False
        if intersecter(regions, r.chrB, r.posB, r.posB + 1):
            br = True
        kind = "extra-regional"
        switch = False
        if not ar and not br:
            if r.chrA == r.chrB and r.posA > r.posB:  # Put non-region first
                switch = True
        if (br and not ar) or (not br and ar):
            kind = "hemi-regional"
            if not br and ar:
                switch = True
        if ar and br:
            if r.chrA == r.chrB:
                rA = regions[r.chrA].overlappingInterval(r.posA, r.posA + 1)
                rB = regions[r.chrB].overlappingInterval(r.posB, r.posB + 1)
                if rA and rB and rA[0] == rB[0] and rA[1] == rB[1]:
                    kind = "intra_regional"
                    if r.posA > r.posB:
                        switch = True
                else:
                    kind = "inter-regional"
                    if r.chrA != sorted([r.chrA, r.chrB])[0]:
                        switch = True
            else:
                kind = "inter-regional"
        if switch:
            chrA, posA, cipos95A, contig2 = r.chrA, r.posA, r.cipos95A, r.contig2
            r.chrA = r.chrB
            r.posA = r.posB
            r.cipos95A = r.cipos95B
            r.chrB = chrA
            r.posB = posA
            r.cipos95B = cipos95A
            r.contig2 = r.contig
            r.contig = contig2
        max_depth = 0
        if kind == "hemi-regional":
            chrom_i = r.chrA
            if chrom_i in regions_depth.chrom_cov_arrays:
                reads_10kb, max_depth = calculate_coverage(r.posA - 10000, r.posA + 10000, regions_depth.chrom_cov_arrays[chrom_i])
                reads_10kb = round(reads_10kb, 3)
            else:
                reads_10kb = 0
        else:
            chrom_i = r.chrA
            if chrom_i in regions_depth.chrom_cov_arrays:
                reads_10kb_left, max_depth = calculate_coverage(r.posA - 10000, r.posA + 10000, regions_depth.chrom_cov_arrays[chrom_i])
                reads_10kb_left = round(reads_10kb_left, 3)
            else:
                reads_10kb_left = 0
            chrom_i = r.chrB
            if chrom_i in regions_depth.chrom_cov_arrays:
                reads_10kb_right, max_depth = calculate_coverage(r.posB - 10000, r.posB + 10000, regions_depth.chrom_cov_arrays[chrom_i])
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
        new_events.append(r)
    return new_events
