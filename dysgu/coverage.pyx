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


def median(L):
    # Stolen from https://github.com/arq5x/lumpy-sv/blob/master/scripts/pairend_distro.py
    if len(L) % 2 == 1:
        return L[int(len(L)/2)]  # cast to int since divisions always return floats in python3
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
    mean = np.median(L) #s / float(len(L))
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
    click.echo("Removed {} outliers with isize >= {}".format(removed, upper_cutoff), err=True)
    mean, stdev = mean_std(L)
    return mean, stdev


class GenomeScanner:

    def __init__(self, inputbam, int max_cov, include_regions, read_threads, buffer_size, regions_only, stdin):

        self.input_bam = inputbam
        self.max_cov = max_cov
        self.include_regions = include_regions
        self.overlap_regions = io_funcs.overlap_regions(include_regions, int_chroms=True, infile=inputbam)
        self.regions_only = regions_only

        self.staged_reads = deque([])
        self.procs = read_threads
        self.buff_size = buffer_size

        self.bam_iter = self.input_bam #inputbam.fetch(until_eof=True)

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

        self.no_tell = True if stdin else False


    def _get_reads(self):

        if not self.include_regions or not self.regions_only:

            while len(self.staged_reads) > 0:
                yield self.staged_reads.popleft()

            tell = 0 if self.no_tell else self.input_bam.tell()
            for a in self.input_bam: # self.bam_iter:

                self._add_to_bin_buffer(a, tell)
                tell = 0 if self.no_tell else self.input_bam.tell()

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

        # Borrowed from lumpy
        cdef int required = 97
        restricted = 3484
        cdef int flag_mask = required | restricted

        cdef int c = 0
        tell = 0
        cdef int flag, tlen
        cdef float approx_read_length


        for a in self.bam_iter:

            tell = 0 if self.no_tell else self.input_bam.tell()
            if self.no_tell:
                self._add_to_bin_buffer(a, tell)
            if len(approx_read_length_l) < 100000:
                flag = a.flag
                if a.seq is not None:
                    rl = a.infer_read_length()
                    if rl:
                        approx_read_length_l.append(rl)
                        if a.rname == a.rnext and flag & flag_mask == required and a.tlen >= 0:

                            inserts.append(a.tlen)

            else:
                break

            if c > 10000000:
                raise ValueError("Cant infer read length after 10 million reads, is max-tlen < 8000?")
            c += 1

        if len(approx_read_length_l) == 0:
            raise RuntimeError("Cant infer read length, no reads?")

        approx_read_length = int(np.median(approx_read_length_l))
        self.approx_read_length = approx_read_length
        if len(inserts) <= 1000 and insert_median == -1:
            insert_median = 300
            insert_stdev = 150

        if insert_median == -1:
            insert_median, insert_stdev = get_insert_params(inserts)

        click.echo(f"Inferred read length {approx_read_length}, "
                   f"insert median {insert_median}, "
                   f"insert stdev {insert_stdev}", err=True)
        self.last_tell = tell

        if not self.no_tell:
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
            click.echo("coverage.pyx No reads, finishing", err=True)
            quit()

        click.echo(f"Total input reads {total_reads}", err=True)

    def add_to_buffer(self, r, n1):

        if self.first == 1:
            # Not possible to get tell position from first read in region, so put into buffer instead
            # if self.procs == 1:
            self.read_buffer[n1] = r
            # else:
            #     self.read_buffer[n1] = r.to_string()

            self.first = 0

        elif len(self.read_buffer) < self.buff_size:
            # if self.procs == 1:
            self.read_buffer[n1] = r
            # else:
            #     self.read_buffer[n1] = r.to_string()
        elif self.no_tell:
            raise BufferError("Read buffer has overflowed, increase --buffer-size")


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

            # click.echo(f"{self.input_bam.get_reference_name(rname)}, coverage bins={ref_length}", err=True)

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
