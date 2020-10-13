#cython: language_level=3, boundscheck=False

from __future__ import absolute_import
from collections import deque, defaultdict
import click
from dysgu import io_funcs
import numpy as np
cimport numpy as np

# from dysgu.map_set_utils cimport robin_set
# from dysgu.map_set_utils cimport hash as xxhash

# from pysam.libchtslib cimport bam_hdr_t, BAM_CINS, BAM_CDEL, BAM_CSOFT_CLIP, \
#     bam_cigar_op, bam_cigar_oplen, bam_get_qname, bam_get_cigar, bam_init1, \
#     htsFile, bam_destroy1, bam_aux_get, sam_read1, bam1_t, bam_dup1, hts_set_threads
# from pysam.libcalignedsegment cimport makeAlignedSegment, AlignedSegment
# from pysam.libcalignmentfile cimport AlignmentFile, AlignmentHeader

# from cython.operator import dereference

# from libc.stdint cimport uint16_t, uint32_t, uint64_t
# from libcpp.deque cimport deque as cpp_deque
# from libcpp.pair cimport pair as cpp_pair

# ctypedef cpp_pair[uint64_t, bam1_t*] deque_item
DTYPE = np.float
ctypedef np.float_t DTYPE_t


def echo(*args):
    click.echo(args, err=True)


def merge_intervals(intervals, srt=True, pad=0, add_indexes=False):
    """
    >>> merge_intervals( [('chr1', 1, 4), ('chr1', 2, 5), ('chr2', 3, 5)] )
    >>> [['chr1', 1, 5], ['chr2', 3, 5]]
    """
    if srt:
        sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))  # by chrom, start, end (index)
    else:
        sorted_by_lower_bound = intervals

    if pad:
        if not add_indexes:
            sorted_by_lower_bound = [[c, 0 if i - pad < 0 else i - pad, j + pad] for c, i, j in sorted_by_lower_bound]
        else:
            sorted_by_lower_bound = [[c, 0 if i - pad < 0 else i - pad, j + pad, k] for c, i, j, k in sorted_by_lower_bound]

    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            if not add_indexes:
                merged.append(higher)
            else:
                merged.append(list(higher)[:3] + [[higher[3]]])
            continue
        elif higher[0] != merged[-1][0]:  # Dont merge intervals on different chroms
            if not add_indexes:
                merged.append(higher)
            else:
                merged.append(list(higher)[:3] + [[higher[3]]])
        else:
            lower = merged[-1]  # Last item on merged (end of interval)
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[1] <= lower[2]:
                if not add_indexes:
                    merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]))
                else:
                    merged[-1] = (lower[0], lower[1], max(higher[2], lower[2]), lower[3] + [higher[3]])
            else:
                if not add_indexes:
                    merged.append(higher)
                else:
                    merged.append(list(higher)[:3] + [[higher[3]]])
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
    click.echo("Removed {} outliers with isize >= {}".format(removed, upper_cutoff), err=True)
    mean, stdev = mean_std(L)
    return mean, stdev


# cdef struct DequeItem:
#     uint64_t hash_val
#     bam1_t* bam_ptr
#     long tell
#
#
# cdef DequeItem makeItem(uint64_t hash_val, bam1_t* bam_ptr, long tell) nogil:
#     cdef DequeItem v
#     v.hash_val = hash_val
#     v.bam_ptr = bam_ptr
#     v.tell = tell
#     return v
#
#
# cdef AlignedSegment passAlignedSegment(bam1_t *src,
#                                        AlignmentHeader header):
#     '''return an AlignedSegment object constructed from `src`'''
#     # note that the following does not call __init__
#     cdef AlignedSegment dest = AlignedSegment.__new__(AlignedSegment)
#     dest._delegate = src #bam_dup1(src)
#     dest.header = header
#     return dest
#
#
# def htslib_iterator(AlignmentFile infile, uint32_t min_within_size, uint32_t clip_length):
#     # Similar to sv2bam.pyx script, but yields SV alignments and tell position rather than outputting to file
#
#     cdef int result
#     cdef htsFile *fp_in = infile.htsfile #hts_open(infile, "r")
#
#     result = hts_set_threads(fp_in, 4)
#     if result != 0:
#         raise IOError("Failed to set threads")
#
#     cdef AlignmentHeader pysam_header = infile.header
#     cdef bam_hdr_t* samHdr = infile.header.ptr  #sam_hdr_read(fp_in)  # read header
#     cdef bam1_t* aln = bam_init1()  # initialize an alignment
#     if not samHdr:
#         raise IOError("Failed to read input header")
#
#     cdef uint32_t max_scope = 100000
#     cdef uint64_t total = 0
#
#     # cdef cpp_pair[uint64_t, bam1_t*] scope_item
#     cdef cpp_deque[DequeItem] scope
#     cdef robin_set[uint64_t] read_names
#
#     cdef uint16_t flag
#     cdef uint32_t* cigar
#     cdef uint32_t k, op, length
#     cdef uint64_t precalculated_hash
#     cdef bam1_t* bam_ptr
#     cdef bam1_t bamt
#
#     cdef long tell = infile.tell()
#     cdef DequeItem scope_item
#     yield_queue = []
#     while sam_read1(fp_in, samHdr, aln) >= 0:
#
#         if scope.size() > max_scope:
#             scope_item = scope[0]
#
#             if read_names.find(scope_item.hash_val, scope_item.hash_val) != read_names.end():
#                 # if not yielded:
#                 yield passAlignedSegment(scope_item.bam_ptr, pysam_header), scope_item.tell
#                     # yielded = 1
#                 # else:
#                 #     bam_destroy1(scope_item.bam_ptr)
#             else:
#                 bam_destroy1(scope_item.bam_ptr)
#             # bam_destroy1(scope_item.bam_ptr)
#             # free(scope_item)
#             scope.pop_front()
#
#         # Process current alignment
#         flag = dereference(aln).core.flag
#         if flag & 1284 or dereference(aln).core.n_cigar == 0 or dereference(aln).core.l_qname == 0:
#             continue
#
#         precalculated_hash = xxhash(bam_get_qname(aln), dereference(aln).core.l_qname, 0)
#         scope.push_back(makeItem(precalculated_hash, bam_dup1(aln), tell))
#         tell = infile.tell() #hts_utell(fp_in)
#         # tell = 1
#         if read_names.find(precalculated_hash, precalculated_hash) == read_names.end():
#             # Check for discordant of supplementary
#             if ~flag & 2 or flag & 2048:
#                 read_names.insert(precalculated_hash)
#                 continue
#
#             # Check for SA tag
#             if bam_aux_get(aln, "SA"):
#                 read_names.insert(precalculated_hash)
#                 continue
#
#             cigar = bam_get_cigar(aln)
#             for k in range(0, dereference(aln).core.n_cigar):
#
#                 op = bam_cigar_op(cigar[k])
#                 if op == -1:
#                     break
#                 length = bam_cigar_oplen(cigar[k])
#                 if length == -1:
#                     break
#
#                 if (op == BAM_CSOFT_CLIP ) and (length >= clip_length):  # || op == BAM_CHARD_CLIP
#                     read_names.insert(precalculated_hash)
#                     break
#
#                 if (op == BAM_CINS or op == BAM_CDEL) and (length >= min_within_size):
#                     read_names.insert(precalculated_hash)
#                     break
#
#     while scope.size() > 0:
#         scope_item = scope[0]
#         if read_names.find(scope_item.hash_val, scope_item.hash_val) != read_names.end():
#             yield passAlignedSegment(scope_item.bam_ptr, pysam_header), scope_item.tell
#         else:
#             bam_destroy1(scope_item.bam_ptr)
#         scope.pop_front()
#
#     # result = hts_close(fp_in)
#     # if result != 0:
#     #     raise IOError("Problem closing input file handle")
#
#     bam_destroy1(aln)


class GenomeScanner:

    def __init__(self, inputbam, int max_cov, include_regions, read_threads, buffer_size, regions_only, stdin,
                 clip_length=30, min_within_size=30):

        self.input_bam = inputbam
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

        self.first = 1  # Not possible to get first position in a target fetch region, so buffer first read instead
        self.read_buffer = dict()

        self.approx_read_length = -1
        self.last_tell = 0
        self.no_tell = True if stdin else False
        self.extended_tags = False

    def _get_reads(self):
        # Two options, reads are collected from whole genome, or from target regions only

        # Scan whole genome
        if not self.include_regions or not self.regions_only:

            # Some reads may have been staged from getting read_length (if file is streamed from stdin)
            while len(self.staged_reads) > 0:
                yield self.staged_reads.popleft()

            tell = 0 if self.no_tell else self.input_bam.tell()
            for aln in self.input_bam:

                self._add_to_bin_buffer(aln, tell)
                tell = 0 if self.no_tell else self.input_bam.tell()

                while len(self.staged_reads) > 0:
                    yield self.staged_reads.popleft()

            if len(self.current_bin) > 0:
                yield self.current_bin

        # Scan input regions
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

            seen_reads = set([])  # Avoid reading same alignment twice, shouldn't happen anyway
            click.echo("Collecting mate-pairs +/- 1kb", err=True)
            for c, s, e in itv:

                tell = -1
                for a in self.input_bam.fetch(c, int(s), int(e)):

                    name = a.qname.__hash__(), a.flag, a.pos
                    if name in seen_reads:
                        continue

                    self._add_to_bin_buffer(a, tell)
                    seen_reads.add(name)

                    while len(self.staged_reads) > 0:
                        yield self.staged_reads.popleft()

                if len(self.current_bin) > 0:
                    yield self.current_bin

    def get_read_length(self, int max_tlen, int insert_median, int insert_stdev, int read_len, ibam=None):
        # This is invoked first to scan the first part of the file for the insert size metrics,
        # or open and process the --ibam alignment file
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

        # Check if tags from dodi have been added to input reads --> add to vcf output if True
        tags_checked = False

        # Use ibam for insert params if possible
        if ibam is not None:
            file_iter = ibam
        else:
            file_iter = self.input_bam

        for a in file_iter:  # self.input_bam:

            if ibam is None:
                tell = 0 if self.no_tell else self.input_bam.tell()
                if self.no_tell:  # Buffer reads if coming from stdin
                    self._add_to_bin_buffer(a, tell)

            if len(approx_read_length_l) < 200000:
                flag = a.flag
                if a.seq is not None:
                    rl = a.infer_read_length()
                    if rl:
                        approx_read_length_l.append(rl)
                        if a.rname == a.rnext and flag & flag_mask == required and a.tlen >= 0:
                            inserts.append(a.tlen)
                            if not tags_checked:
                                if a.has_tag("DP"):
                                    self.extended_tags = True
                                tags_checked = True

            else:
                break

            if c > 20000000:
                raise ValueError("Cant infer read length after 10 million reads, is max-tlen < 8000?")
            c += 1

        if len(approx_read_length_l) == 0:
            raise RuntimeError("Cant infer read length, no reads?")

        approx_read_length = int(np.median(approx_read_length_l))
        self.approx_read_length = approx_read_length
        if len(inserts) <= 100 and insert_median == -1:
            insert_median = 300
            insert_stdev = 150

        if insert_median == -1:
            insert_median, insert_stdev = get_insert_params(inserts)

        click.echo(f"Inferred read length {approx_read_length}, "
                   f"insert median {insert_median}, "
                   f"insert stdev {insert_stdev}", err=True)

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

        # Add last seen coverage bin
        if total_reads == 0:
            click.echo("coverage.pyx found no reads, finishing", err=True)
            quit()

        click.echo(f"Total input reads {total_reads}", err=True)

    def add_to_buffer(self, r, n1):

        if self.first == 1:
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

        # PCR duplicate, fails quality check, read unmapped
        if flag & 1540 or a.cigartuples is None or a.seq is None:
            return

        cdef int rname = a.rname
        cdef int apos = a.pos
        cdef int bin_pos = int(apos / 100)
        cdef int ref_length
        cdef str reference_name = ""
        cdef int aend = a.reference_end
        cdef float current_coverage

        if rname not in self.depth_d:

            # Get the chromosome size from infile
            ref_length = int(self.input_bam.get_reference_length(
                self.input_bam.get_reference_name(rname)) / 100)

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


cdef float add_coverage(int start, int end, DTYPE_t[:] chrom_depth) nogil:

    # Round start and end to get index
    cdef float fs = start / 100
    cdef float fe = end / 100
    cdef int bin_start = <int> fs  # Cast to int
    cdef int bin_end = <int> fe
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


cpdef calculate_coverage(int start, int end, DTYPE_t[:] chrom_depth):
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

#
# support_bands = [0.9, 0.25]
#
# cpdef float adaptive_support_threshold(bunch, depth_table, allele_fraction):
#
#     pad = 500
#     max_val = -1
#
#     positions = []
#
#     bin_counts = defaultdict(int)
#     background_counts = {}
#     for r_info in bunch:
#
#         event_pos = r_info[6]
#
#         # Get bin value
#         p = int(event_pos / 100)
#         bin_counts[p] += 1
#
#         if p not in background_counts:
#             counts, _ = calculate_coverage(event_pos, event_pos, depth_table[r_info[3]])
#             background_counts[p] = counts
#         positions.append(event_pos)
#
#     max_support = len(bunch)
#     if max_support > len(support_bands) - 1:
#         thresh = allele_fraction
#     else:
#         thresh = support_bands[max_support]
#     for k, v in background_counts.items():
#         if v == 0:
#             continue
#         af = float(bin_counts[k]) / v
#         if af > max_val:
#             max_val = af
#         if af > thresh:
#             return 1
#
#     if max_val == -1:
#         return 1
#     else:
#         return 0


cpdef dict get_raw_coverage_information(r, regions, regions_depth, infile, max_cov):

    # Check if side A in regions
    ar = False  # c_io_funcs.intersecter_int_chrom
    if io_funcs.intersecterpy(regions, r["chrA"], r["posA"], r["posA"] + 1):
        ar = True
    br = False
    if io_funcs.intersecterpy(regions, r["chrB"], r["posB"], r["posB"] + 1):
        br = True

    kind = "extra-regional"
    if not ar and not br:
        if r["chrA"] == r["chrB"] and r["posA"] > r["posB"]:  # Put non-region first
            switch = True

        # Skip if regions have been provided; almost always false positives?
        # if regions is not None:
        #     return None

    switch = False
    if (br and not ar) or (not br and ar):
        kind = "hemi-regional"
        if not br and ar:
            switch = True

    if ar and br:
        if r["chrA"] == r["chrB"]:
            rA = list(regions[r["chrA"]].find_overlap(r["posA"], r["posA"] + 1))[0]
            rB = list(regions[r["chrB"]].find_overlap(r["posB"], r["posB"] + 1))[0]
            if rA[0] == rB[0] and rA[1] == rB[1]:
                kind = "intra_regional"
                # Put posA first
                if r["posA"] > r["posB"]:
                    switch = True
            else:
                kind = "inter-regional"
                if r["chrA"] != sorted([r["chrA"], r["chrB"]])[0]:
                    switch = True
        else:
            kind = "inter-regional"

    if switch:
        chrA, posA, cipos95A, contig2 = r["chrA"], r["posA"], r["cipos95A"], r["contig2"]
        r["chrA"] = r["chrB"]
        r["posA"] = r["posB"]
        r["cipos95A"] = r["cipos95B"]
        r["chrB"] = chrA
        r["posB"] = posA
        r["cipos95B"] = cipos95A
        r["contig2"] = r["contig"]
        r["contig"] = contig2

    max_depth = 0
    if kind == "hemi-regional":
        chrom_i = infile.get_tid(r["chrA"])
        if chrom_i in regions_depth:
            reads_10kb, max_depth = calculate_coverage(r["posA"] - 10000, r["posA"] + 10000, regions_depth[chrom_i])
            reads_10kb = round(reads_10kb, 3)
        else:
            reads_10kb = 0
    else:
        # Calculate max
        chrom_i = infile.get_tid(r["chrA"])
        if chrom_i in regions_depth:
            reads_10kb_left, max_depth = calculate_coverage(r["posA"] - 10000, r["posA"] + 10000, regions_depth[chrom_i])
            reads_10kb_left = round(reads_10kb_left, 3)
        else:
            reads_10kb_left = 0
        chrom_i = infile.get_tid(r["chrB"])
        if chrom_i in regions_depth:
            reads_10kb_right, max_depth = calculate_coverage(r["posB"] - 10000, r["posB"] + 10000, regions_depth[chrom_i])
            reads_10kb_right = round(reads_10kb_right, 3)
        else:
            reads_10kb_right = 0
        if reads_10kb_left > reads_10kb_right:
            reads_10kb = reads_10kb_left
        else:
            reads_10kb = reads_10kb_right

    if reads_10kb >= max_cov * 0.55 and not ar and not br:
        return None

    r["kind"] = kind
    r["raw_reads_10kb"] = reads_10kb
    if r["chrA"] != r["chrB"]:
        r["svlen"] = 1000000

    r["mcov"] = max_depth
    return r
