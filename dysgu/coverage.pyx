#cython: language_level=3, boundscheck=False

from __future__ import absolute_import
from collections import deque
import click
from dysgu import io_funcs
import numpy as np
cimport numpy as np
import mmh3

import pysam
# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
from libc.stdint cimport uint32_t, uint8_t


DTYPE = np.float
ctypedef np.float_t DTYPE_t


def echo(*args):
    click.echo(args, err=True)


# class PyAlignment(object):
#     """Picklable struct to hold the contents of pysam alignment"""
#     __slots__ = ["reference_end", "cigar", "pos", "flag", "rname", "qname", "rnext", "pnext", "seq",
#                  "mapq", "cigartuples",
#                  "query_qualities", "has_SA", "tlen", "has_NM", "has_DP", "has_DN", "has_DA", "has_NP"]
#
#     def __init__(self, a):
#         self.reference_end = a.reference_end
#         self.cigar = a.cigar
#         self.pos = a.pos
#         self.flag = a.flag
#         self.rname = a.rname
#         self.qname = a.qname
#         self.rnext = a.rnext
#         self.pnext = a.pnext
#         self.seq = a.seq
#         self.mapq = a.mapq
#         self.cigartuples = self.cigar
#         self.query_qualities = a.query_qualities
#         self.has_SA = a.has_tag("SA")
#         self.tlen = a.tlen
#         self.has_NM = False
#         self.has_DP = False
#         self.has_DN = False
#         self.has_DA = False
#         self.has_NP = False
#
#         for key in ("NM", "DP", "DN", "DA", "NP"):
#             has = a.has_tag(key)
#             if has:
#                 setattr(self, f"has_{key}", a.get_tag(key))
#
#     def has_tag(self, tag):
#         if getattr(self, f"has_{tag}", None):
#             return True
#         else:
#             return False
#     def get_tag(self, tag):
#         return getattr(self, f"has_{tag}", None)


def merge_intervals(intervals, srt=True, pad=0):
    """
    >>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [['chr1', 1, 5], ['chr2', 3, 5]]
    """
    if srt:
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))  # by chrom, start, end
    else:
        sorted_by_lower_bound = intervals

    if pad:
        sorted_by_lower_bound = [(c, 0 if i - pad < 0 else i - pad, j + pad) for c, i, j in intervals]

    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
            continue
        elif higher[0] != merged[-1][0]:  # Dont merge intervals on different chroms
            merged.append(higher)
        else:
            lower = merged[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= lower[2]:
                merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]))
            else:
                merged.append(higher)
    return merged


class GenomeScanner:

    def __init__(self, inputbam, int max_cov, include_regions, read_threads, buffer_size, regions_only):

        self.input_bam = inputbam
        self.max_cov = max_cov
        self.include_regions = include_regions
        self.overlap_regions = io_funcs.overlap_regions(include_regions, int_chroms=True, infile=inputbam)
        self.regions_only = regions_only

        self.staged_reads = deque([])
        self.procs = read_threads
        self.buff_size = buffer_size

        self.bam_iter = inputbam.fetch(until_eof=True)

        self.current_bin = []
        self.current_cov = 0
        self.current_chrom = 0
        self.current_pos = 0
        self.current_cov_array = None

        self.reads_dropped = 0
        self.depth_d = {} # Use array for fast lookup # defaultdict(lambda: defaultdict(int))

        self.first = 1  # Not possible to get first position in a target fetch region, so buffer first read instead
        self.read_buffer = dict()

        self.approx_read_length = -1

        self.last_tell = 0

    def _get_reads(self):

        if not self.include_regions or not self.regions_only:

            tell = self.input_bam.tell()
            for a in self.bam_iter:

                self._add_to_bin_buffer(a, tell)
                tell = self.input_bam.tell()

                while len(self.staged_reads) > 0:
                    yield self.staged_reads.popleft()

            if len(self.current_bin) > 0:
                yield self.current_bin

        else:
            # Reads must be fed into graph in sorted order, find regions of interest first
            intervals_to_check = []  # Containing include_regions, and also mate pairs
            pad = 1000

            regions = [i.strip().split("\t")[:3] for i in open(self.include_regions, "r") if i[0] != "#"]
            for c, s, e in regions:
                intervals_to_check.append((c, int(s), int(e)))

            for c, s, e in regions:

                for a in self.input_bam.fetch(c, int(s), int(e)):
                    # Mate unmapped, not primary, fails QC, duplicate
                    if not a.flag & 1800:
                        p1 = a.pnext - pad
                        c2 = self.input_bam.get_reference_name(a.rnext)
                        intervals_to_check.append((c2, 0 if p1 < 0 else p1, a.pnext + pad))
                        if a.has_tag("SA"):
                            sa = a.get_tag("SA").split(";")
                            for v in sa:
                                chrom2, pos2, _ = v.split(",", 2)
                                pos2 = int(pos2) - pad
                                intervals_to_check.append((chrom2, 0 if pos2 < 0 else pos2, pos2 + pad))

            itv = merge_intervals(intervals_to_check)

            seen_reads = set([])  # Avoid reading same alignment twice
            click.echo("Collecting mate-pairs +/- 1kb", err=True)
            for c, s, e in itv:

                tell = -1
                for a in self.input_bam.fetch(c, int(s), int(e)):

                    name = mmh3.hash(a.qname, 42), a.flag, a.pos
                    if name in seen_reads:
                        continue

                    self._add_to_bin_buffer(a, tell)
                    seen_reads.add(name)

                    while len(self.staged_reads) > 0:
                        yield self.staged_reads.popleft()

                if len(self.current_bin) > 0:
                    yield self.current_bin

    def get_read_length(self, int max_tlen, int insert_median, int insert_stdev, int read_len):

        if read_len != -1:
            click.echo(f"Read length {read_len}, "
                   f"insert_median {insert_median}, "
                   f"insert stdev {insert_stdev}", err=True)
            self.approx_read_length = read_len
            return insert_median, insert_stdev

        approx_read_length_l = []
        inserts = []

        cdef int c = 0
        tell = 0
        cdef int flag, tlen
        cdef float approx_read_length

        for a in self.bam_iter:

            tell = self.input_bam.tell()
            if len(approx_read_length_l) < 20000:
                flag = a.flag
                if a.seq is not None:

                    if not flag & 3852:
                        approx_read_length_l.append(a.infer_read_length())
                        tlen = abs(a.tlen)
                        if not flag & 2 and 0 < tlen < max_tlen:  # Skip discordants  # Todo donsnt work for single end
                            if insert_median == -1:
                                inserts.append(tlen)
            else:
                break

            if c > 10000000:
                raise ValueError("Cant infer read length after 10 million reads, is max-tlen < 8000?")
            c += 1

        if len(approx_read_length_l) == 0:
            raise RuntimeError("Cant infer read length, no reads?")

        approx_read_length = int(np.median(approx_read_length_l))
        self.approx_read_length = approx_read_length
        if len(inserts) == 0 and insert_median == -1:
            insert_median = 300
            insert_stdev = 150

        if insert_median == -1:
            insert_stdev = int(np.std(inserts))
            insert_median = int(np.median(inserts))

        click.echo(f"Inferred Read length {approx_read_length}, "
                   f"insert_median {insert_median}, "
                   f"insert stdev {insert_stdev}", err=True)
        self.last_tell = tell

        self.input_bam.reset()
        return insert_median, insert_stdev

    def iter_genome(self):
        # Send any staged reads
        cdef int total_reads = 0

        for staged in self._get_reads():
            total_reads += len(staged)
            yield staged

        # Add last seen coverage bin
        if total_reads == 0:
            click.echo("No reads, finishing", err=True)
            quit()

        click.echo(f"Total input reads {total_reads}", err=True)

    def add_to_buffer(self, r, n1):

        if self.first == 1:
            # Not possible to get tell position from first read in region, so put into buffer instead
            if self.procs == 1:
                self.read_buffer[n1] = r
            else:
                self.read_buffer[n1] = r.to_string()

            self.first = 0

        elif len(self.read_buffer) < self.buff_size:
            if self.procs == 1:
                self.read_buffer[n1] = r
            else:
                self.read_buffer[n1] = r.to_string()


    def _add_to_bin_buffer(self, a, tell):

        cdef int flag = a.flag

        # PCR duplicate, fails quality check, not primary alignment, read unmapped
        if flag & 1796 or a.cigartuples is None or a.seq is None:
            return

        cdef int rname = a.rname
        cdef int apos = a.pos
        cdef int bin_pos = int(apos / 100)
        cdef int ref_length
        cdef str reference_name = ""
        cdef int aend = a.infer_query_length() + apos
        cdef float current_coverage

        if rname not in self.depth_d:

            # Get the chromosome size from infile
            ref_length = int(self.input_bam.get_reference_length(
                self.input_bam.get_reference_name(rname)) / 100)

            click.echo(f"{self.input_bam.get_reference_name(rname)}, coverage bins={ref_length}", err=True)

            # Define a big numpy array to hold count information
            self.depth_d[rname] = np.zeros(ref_length + 1, dtype=np.float)

        if self.current_chrom != rname:
            self.current_chrom = rname
            self.current_cov_array = self.depth_d[rname]
        elif self.current_cov_array is None:
            self.current_cov_array = self.depth_d[rname]

        current_cov = add_coverage(apos, aend, self.current_cov_array)

        in_roi = False

        if self.overlap_regions:
            in_roi = io_funcs.intersecter_int_chrom(self.overlap_regions, a.rname, apos, apos+1)

        if rname == self.current_chrom and bin_pos == self.current_pos:

            if current_cov >= self.max_cov and not in_roi:
                if len(self.current_bin) > 0:
                    self.current_bin = []
                    self.reads_dropped += len(self.current_bin)
                self.reads_dropped += 1
                return

            self.current_bin.append((a, tell))  # Add to staging area

        else:  # New staged bin

            if len(self.current_bin) != 0 and (current_cov < self.max_cov or in_roi):
                self.staged_reads.append(self.current_bin)  # Send for further processing

            self.current_chrom = rname
            self.current_pos = bin_pos
            self.current_bin = [(a, tell)]


cpdef float add_coverage(int start, int end, DTYPE_t[:] chrom_depth) nogil:

    # Round start and end to get index
    cdef float fs = start / 100
    cdef float fe = end / 100

    cdef int bin_start = <int> fs  # Cast to int
    cdef int bin_end = <int > fe

    if bin_start < 0:
        bin_start = 0
    if bin_end > <int> len(chrom_depth) - 1:
        bin_end = <int> len(chrom_depth) - 1

    # Fraction overlapping start and end bins
    cdef float ol_start = <float> bin_start + 1 - fs

    chrom_depth[bin_start] += ol_start

    cdef float ol_end = 0
    if bin_start != bin_end:  # Read spans more than one bin
        ol_end = fe - (<float> bin_end)

        chrom_depth[bin_end] += ol_end

        if bin_end - bin_start > 1:
            # Fill in between
            for i in range(bin_start + 1, bin_end):
                chrom_depth[i] += 1

    return chrom_depth[bin_start]


# def scan_whole_genome(inputbam, int max_cov, tree):
#
#     chrom_lengths = {inputbam.get_reference_name(v): k for v, k in enumerate(inputbam.lengths)}
#     chrom_order = deque([inputbam.get_reference_name(v) for v, k in enumerate(inputbam.lengths)])
#
#     depth_d = defaultdict(lambda: defaultdict(int))
#     approx_read_length_l = []
#
#     cdef int bin_pos
#     cdef int bin_reads = 0
#     cdef int rname
#     cdef float approx_read_length = -1
#     cdef int total_dropped = 0
#     cdef int reads_dropped = 0
#     cdef int current_chrom = 0
#     cdef int current_pos = 0
#     cdef float current_cov = 0.
#
#     bad_intervals = []
#
#     c = 0
#     for a in inputbam.fetch(until_eof=True):
#
#         # fails quality checks, PCR duplicate, unmapped, mate unmapped, not primary, proper pair, supplementary
#         # read unmapped, not primary, fails checks, optical duplicate, supplementary
#         if a.flag & 3844:
#             continue
#
#         if approx_read_length == -1:
#             if len(approx_read_length_l) < 10000:
#                 if a.seq is not None and not a.flag & 2048:
#                     approx_read_length_l.append(a.infer_read_length())
#             else:
#                 approx_read_length = np.median(approx_read_length_l)
#                 assert approx_read_length != -1
#
#         rname = a.rname
#         bin_pos = int((int(a.pos) / 100)) * 100
#
#         if rname == current_chrom and bin_pos == current_pos:
#
#             if current_cov > max_cov:
#                 reads_dropped += 1
#                 if approx_read_length == -1:
#                     current_cov += (a.infer_read_length() / 100.)
#                 else:
#                     current_cov += (approx_read_length / 100.)
#
#                 continue  # Already added to bad intervals below
#
#             if approx_read_length == -1:
#                 current_cov += (a.infer_read_length() / 100.)
#             else:
#                 current_cov += (approx_read_length / 100.)
#
#
#         else:
#             depth_d[inputbam.get_reference_name(current_chrom)][current_pos] = current_cov
#
#             current_chrom = rname
#             current_pos = bin_pos
#             if approx_read_length == -1:
#                 current_cov = (a.infer_read_length() / 100.)
#             else:
#                 current_cov = (approx_read_length / 100.)
#             bin_reads = 1
#
#         if current_cov > max_cov:
#             bad_intervals.append((inputbam.get_reference_name(current_chrom), current_pos, current_pos + 101))
#             reads_dropped += bin_reads
#             total_dropped += 100
#
#
#     # Add last bin
#     if current_cov < max_cov:
#         depth_d[inputbam.get_reference_name(current_chrom)][current_pos] = current_cov
#
#     approx_read_length = np.median(approx_read_length_l)
#     merged_bad = merge_intervals(bad_intervals, srt=False, pad=1)
#
#     complemented_intervals = complement(chrom_lengths, chrom_order, merged_bad)
#     click.echo("Skipping {} kb of reference, primary read-coverage >= {}X ({} reads)".format(
#         round(total_dropped / 1e3, 3), max_cov, reads_dropped), err=True)
#
#     return depth_d, approx_read_length, complemented_intervals
#
#
# def scan_regions(inputbam, include_, max_cov, tree):
#     # Scan the input regions, then re-scan the locations of the mate-pairs around the genome
#     depth_d = {}
#
#     roi = data_io.get_include_reads(include_, inputbam)
#     approx_read_length = []
#     inserts = []
#
#     intervals_to_check = set([])  # A list of regions containing many bins
#     target_regions = set([])  # A list of individual bins to work on
#     for a in roi:
#
#         if a.flag & 1540:
#             # Unmapped, fails Q, duplicate
#             continue
#
#         if a.flag & 266:
#             # not primary, mate unmapped, proper-pair
#             continue
#
#         # Skip if both positions are non-canonical chromosomes
#         rname, rnext = inputbam.get_reference_name(a.rname), inputbam.get_reference_name(a.rnext)
#
#         bin_pos1 = int((int(a.pos) / 100)) * 100
#         bin_pos2 = int((int(a.pnext) / 100)) * 100
#
#         # Need to check 10kb around pair and mate, find all intervals to check (limits where SVs will be called)
#         c1 = (rname, 0 if bin_pos1 - 10000 < 0 else bin_pos1 - 10000, bin_pos1 + 10000)
#         c2 = (rnext, 0 if bin_pos2 - 10000 < 0 else bin_pos2 - 10000, bin_pos2 + 10000)
#
#         intervals_to_check.add(c1)
#         intervals_to_check.add(c2)
#
#         if len(approx_read_length) < 1000:
#             if not a.flag & 2048:
#                 approx_read_length.append(len(a.infer_read_length()))
#                 if abs(a.template_length) < 1000:
#                     inserts.append(int(abs(a.template_length)))
#
#         # Target regions include upstream and downstream of target - means supplementary reads are not missed later
#         target_regions.add((rname, bin_pos1))
#         target_regions.add((rnext, bin_pos2))
#
#     # Get coverage within 10kb
#     merged_to_check = merge_intervals(intervals_to_check, srt=True)
#
#     # Get coverage in regions
#     for item in merged_to_check:
#         for a in inputbam.fetch(*item):
#
#             if a.flag & 1540:
#                 # Unmapped, fails Q, duplicate
#                 continue
#
#             rname = item[0]
#             bin_pos1 = int(int(a.pos) / 100) * 100
#
#             # Get depth at this locus
#             # Use an ordered dict, obviates need for sorting when merging intervals later on
#             if rname not in depth_d:
#                 depth_d[rname] = SortedDict()
#             if bin_pos1 not in depth_d[rname]:
#                 depth_d[rname][bin_pos1] = 1
#             else:
#                 depth_d[rname][bin_pos1] += 1
#
#     approx_read_length = np.array(approx_read_length)
#     if len(inserts) > 0:
#         insert_std = np.array(inserts).std()
#     else:
#         insert_std = 600
#     approx_read_length = approx_read_length.mean()
#
#     # Increase target region space here. Currently target regions only point to primary alignments, increase to capture
#     # supplementary mappings nearby
#     pad = insert_std * 3
#
#     new_targets = set([])
#     for chrom, primary_site in target_regions:
#         # Pad is determined by insert stdev
#         lower_bin = int((primary_site - pad) / 100) * 100
#         lower_bin = 0 if lower_bin < 0 else lower_bin
#         upper_bin = (int((primary_site + pad) / 100) * 100) + 100
#         for s in range(lower_bin, upper_bin, 100):
#             new_targets.add((chrom, s))
#
#     total_dropped, reads_dropped = 0, 0
#     intervals = []
#     for chrom, bin_start in new_targets:
#
#         if bin_start in depth_d[chrom]:
#             d = depth_d[chrom][bin_start]
#         else:
#             d = 0  # No reads found
#         bin_cov = (d * approx_read_length) / 100
#
#         if bin_cov < max_cov or data_io.intersecter(tree, chrom, bin_start, bin_start + 100):
#             intervals.append((chrom, int(bin_start), int(bin_start + 101)))  # 1 bp overlap with adjacent window
#         else:
#             total_dropped += 100
#             reads_dropped += d
#
#     click.echo("Skipping {} kb of reference, read-coverage >= {}X ({} reads)".format(
#         round(total_dropped / 1e3, 3), max_cov, reads_dropped), err=True)
#
#     merged_targets = merge_intervals(intervals, srt=True)
#
#     return depth_d, approx_read_length, merged_targets
#
#
# def get_low_coverage_regions(inputbam, max_cov, tree, include_):
#     # depth_d contains read counts +- 10kb of regions of interest, merged is used to construct the main graph
#     if tree is None:
#         depth_d, approx_read_length, merged = scan_whole_genome(inputbam, max_cov, tree)
#
#     else:
#         depth_d, approx_read_length, merged = scan_regions(inputbam, include_, max_cov, tree)
#
#     return merged, approx_read_length, depth_d


cpdef float calculate_coverage(int start, int end, DTYPE_t[:] chrom_depth) nogil:
    # Round start and end to get index
    cdef float fs = start / 100
    cdef float fe = end / 100

    start = <int> fs  # Cast to int
    end = <int > fe
    if start < 0:
        start = 0
    cdef int len_chrom = <int> len(chrom_depth)
    if end > len_chrom:
        end = len_chrom
    cdef int i
    cdef float total = 0

    for i in range(start, end):
        total += chrom_depth[i]

    if total == 0:
        return 0
    return total / (end - start)
