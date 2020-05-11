#cython: language_level=3

from __future__ import absolute_import
import pysam
import time
import datetime
import numpy as np
import click
from subprocess import call
from collections import deque
import resource


from libc.stdint cimport uint32_t, uint16_t, uint64_t
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair

from dysgu import io_funcs
from dysgu cimport map_set_utils
from dysgu.map_set_utils cimport robin_set
from dysgu.coverage import get_insert_params

from dysgu.map_set_utils cimport hash as xxhash

# from pysam.libchtslib cimport bam_hdr_t, BAM_CINS, BAM_CDEL, BAM_CSOFT_CLIP, \
#     bam_cigar_op, bam_cigar_oplen, bam_get_qname, bam_get_cigar, \
#     hts_open, hts_set_threads, sam_hdr_read, bam_init1, \
#     sam_hdr_write, htsFile, bam_destroy1, hts_close, bam_aux_get, sam_write1, sam_read1, bam1_t, bam_dup1

from cython.operator import dereference

# ctypedef cpp_pair[uint64_t, bam1_t*] deque_item


# cdef int search_alignments(char* infile, char* outfile, uint32_t min_within_size, uint32_t clip_length,
#                            int threads) nogil:
#
#     cdef int result
#     cdef htsFile *fp_in = hts_open(infile, "r")
#     if threads:
#         result = hts_set_threads(fp_in, threads)
#         if result != 0:
#             raise IOError("Failed to set threads")
#
#     cdef bam_hdr_t* samHdr = sam_hdr_read(fp_in)  # read header
#     cdef bam1_t* aln = bam_init1()  # initialize an alignment
#     if not samHdr:
#         raise IOError("Failed to read input header")
#
#     cdef htsFile *f_out = hts_open(outfile, "wb0")
#
#     result = sam_hdr_write(f_out, samHdr)
#     if result != 0:
#         raise IOError("Failed to write header to output")
#
#     cdef uint32_t max_scope = 100000
#     cdef uint64_t total = 0
#
#     cdef cpp_pair[uint64_t, bam1_t*] scope_item
#     cdef cpp_deque[deque_item] scope
#     cdef robin_set[uint64_t] read_names
#
#     cdef uint16_t flag
#     cdef uint32_t* cigar
#     cdef uint32_t k, op, length
#     cdef uint64_t precalculated_hash
#     cdef bam1_t* bam_ptr
#     cdef bam1_t bamt
#
#     while sam_read1(fp_in, samHdr, aln) >= 0:
#
#         if scope.size() > max_scope:
#             scope_item = scope[0]
#
#             if read_names.find(scope_item.first, scope_item.first) != read_names.end():
#                 result = sam_write1(f_out, samHdr, scope_item.second)
#                 if result < 0:
#                     raise IOError("Problem writing alignment record")
#                 total += 1
#             # free(dereference(scope_item.second).data)
#             # free(scope_item.second)
#             bam_destroy1(scope_item.second)
#             scope.pop_front()
#
#         # Process current alignment
#         flag = dereference(aln).core.flag
#         if flag & 1284 or dereference(aln).core.n_cigar == 0 or dereference(aln).core.l_qname == 0:
#             continue
#
#         precalculated_hash = xxhash(bam_get_qname(aln), dereference(aln).core.l_qname, 0)
#
#         bam_ptr = bam_dup1(aln)
#         scope.push_back(deque_item(precalculated_hash, bam_ptr))  # bam_dup1(aln)
#
#         if read_names.find(precalculated_hash, precalculated_hash) == read_names.end():
#
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
#
#             # Check cigar
#             # for (uint32_t k = 0; k < aln->core.n_cigar; k++):
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
#         if read_names.find(scope_item.first, scope_item.first) != read_names.end():
#             result = sam_write1(f_out, samHdr, scope_item.second)
#             if result < 0:
#                 raise IOError("Problem writing alignment record")
#             total += 1
#
#         scope.pop_front()
#
#     result = hts_close(fp_in)
#     if result != 0:
#         raise IOError("Problem closing input file handle")
#
#     result = hts_close(f_out)
#     if result != 0:
#         raise IOError("Problem closing output file handle")
#     f_out = NULL
#
#     bam_destroy1(aln)
#     return total


def echo(*args):
    click.echo(args, err=True)

def iter_bam(bam, search, show=True):

    if not search:
        click.echo("Searching input file", err=True)
        for aln in bam.fetch(until_eof=True):  # Also get unmapped reads
            yield aln
    else:
        if show:
            click.echo("Limiting search to {}".format(search), err=True)

        with open(search, "r") as bed:
            for line in bed:
                if line[0] == "#":
                    continue
                chrom, start, end = line.split("\t")[:3]
                for aln in bam.fetch(reference=chrom, start=int(start), end=int(end)):  # Also get unmapped reads

                    yield aln


def config(args):

    kind = args["bam"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "r"}
    if kind not in opts:
        raise ValueError("Input file format not recognized, use .bam,.sam or .cram extension")

    click.echo("Input file: {}".format(args["bam"], kind), err=True)

    try:
        bam = pysam.AlignmentFile(args["bam"], opts[kind], threads=args["procs"])
    except:
        raise RuntimeError("Problem opening bam/sam/cram, check file has .bam/.sam/.cram in file name, and file has a header")

    bam_i = iter_bam(bam, args["search"])



    clip_length = args["clip_length"]

    send_output = None
    v = ""
    if args["output"] != "None":
        if args["output"] == "stdout":
            v = "-"
        else:
            v = args["output"]
        send_output = pysam.AlignmentFile(v, "wb", template=bam)

    out_name = "-" if (args["reads"] == "-" or args["reads"] == "stdout") else args["reads"]

    reads_out = pysam.AlignmentFile(out_name, "wbu", template=bam)

    # outbam = pysam.AlignmentFile(out_name + ".bam", "wbu", template=bam, threads=args["procs"])

    return bam, bam_i, clip_length, send_output, reads_out


cdef tuple get_reads(bam, bam_i, exc_tree, uint32_t clip_length, send_output, outbam):

    # cdef uint32_t clip_length
    cdef int flag
    # cdef str qname
    cdef long qname
    cdef list cigartuples

    # Borrowed from lumpy
    cdef int required = 97
    restricted = 3484
    cdef int flag_mask = required | restricted

    scope = deque([])

    # read_names = set([])
    cdef robin_set[long] read_names

    insert_size = []
    read_length = []

    cdef int max_scope = 100000
    cdef int count = 0
    cdef int nn = 0
    cdef bint paired_end = 0

    for nn, r in enumerate(bam_i):

        if send_output:
            send_output.write(r)

        while len(scope) > max_scope:
            qname, query = scope.popleft()

            if read_names.find(qname, qname) != read_names.end():

            # if query.qname.__hash__() in read_names:
                outbam.write(query)
                count += 1

        if exc_tree:  # Skip exclude regions
            if io_funcs.intersecter_int_chrom(exc_tree, r.rname, r.pos, r.pos + 1):
                continue

        flag = r.flag
        if flag & 1284 or r.cigartuples is None:
            continue  # Read is duplicate or not primary alignment, or unmapped

        if flag & 1 and len(insert_size) < 100000:  # read_paired
            paired_end = 1
            if r.seq is not None:
                if r.rname == r.rnext and r.flag & flag_mask == required and r.tlen >= 0:
                    read_length.append(r.infer_read_length())
                    insert_size.append(r.tlen)


        qname = r.qname.__hash__()
        scope.append((qname, r))

        # if qname not in read_names:
        if read_names.find(qname, qname) == read_names.end():
            if map_set_utils.cigar_clip(r, clip_length):
                # read_names.add(qname)
                read_names.insert(qname)
            elif ~ flag & 2 or flag & 2048:  # Save if read is discordant or supplementary
                # read_names.add(qname)
                read_names.insert(qname)
            elif r.has_tag("SA"):
                # read_names.add(qname)
                read_names.insert(qname)
            elif any((j == 1 or j == 2) and k >= 30 for j, k in r.cigartuples):
                # read_names.add(qname)
                read_names.insert(qname)

    while len(scope) > 0:
        qname, query = scope.popleft()
        if read_names.find(qname, qname) != read_names.end():
        # if query.qname.__hash__() in read_names:
            outbam.write(query)
            count += 1

    outbam.close()
    if send_output:
        send_output.close()

    click.echo("Total input reads in bam {}".format(nn + 1), err=True)
    insert_median, insert_stdev, approx_read_length = 0, 0, 0

    if paired_end:
        if len(read_length) == 0:
            raise ValueError("No paired end reads")
        approx_read_length = int(np.mean(read_length))
        if len(insert_size) == 0:
            insert_median, insert_stdev = 300, 150  # args["insert_median"], args["insert_stdev"]
            click.echo("WARNING: could not infer insert size, no 'normal' pairings found. Using arbitrary values", err=True)
        else:
            insert_median, insert_stdev = get_insert_params(insert_size)
            insert_median, insert_stdev = np.round(insert_median, 2), np.round(insert_stdev, 2)
        click.echo(f"Inferred read length {approx_read_length}, "
                       f"insert median {insert_median}, "
                       f"insert stdev {insert_stdev}", err=True)

    return count, insert_median, insert_stdev, approx_read_length


cdef extern from "find_reads.h":
    cdef int search_hts_alignments(char* infile, char* outfile, uint32_t min_within_size, uint32_t clip_length,
                                    int threads)

def process(args):
    t0 = time.time()

    # cdef int paired_end = int(args["paired"] == "True")

    exc_tree = None
    if args["exclude"]:
        click.echo("Excluding {} from search".format(args["exclude"]), err=True)
        exc_tree = io_funcs.overlap_regions(args["exclude"])

    out_name = "-" if (args["reads"] == "-" or args["reads"] == "stdout") else args["reads"]

    cdef bytes infile_string = args["bam"].encode("ascii")
    cdef bytes outfile_string = out_name.encode("ascii")

    if args["output"] == "None" and exc_tree is None:  # and paired_end:

        t0 = time.time()

        count = search_hts_alignments(infile_string, outfile_string, 30, args["clip_length"], args["procs"])
        if count < 0:
            click.echo("Error reading input file", err=True)
            quit()
        echo(time.time() - t0)
        click.echo("dysgu fetch {}, n={}, mem={} Mb, time={} h:m:s".format(args["bam"],
                                                                    count,
                                                                   int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
                                                                   str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)
        quit()
        return {}

    else:
        echo("pysam version")
        bam, bam_i, clip_length, send_output, outbam = config(args)
        count, insert_median, insert_stdev, read_length = get_reads(bam, bam_i, exc_tree, clip_length, send_output, outbam,
                                                                )


    if args["index"] == "True" and args["reads"] not in "-,stdout":
        call("samtools index -@{} {}".format(args["procs"], args["reads"]), shell=True)

    click.echo("dysgu fetch {}, n={}, mem={} Mb, time={} h:m:s".format(args["bam"],
                                                                    count,
                                                                   int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
                                                                   str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)
    click.echo(time.time() - t0, err=True)
    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            }

