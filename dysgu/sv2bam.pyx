#cython: language_level=3

from __future__ import absolute_import
import pysam
import time
import datetime
import numpy as np
import os
import click
from subprocess import call
from collections import deque
import resource

from libc.stdint cimport uint32_t

from dysgu import data_io, io_funcs
from dysgu cimport map_set_utils
from dysgu.coverage import get_insert_params


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
    opts = {"bam": "rb", "cram": "rc", "sam": "rs"}
    if kind not in opts:
        raise ValueError("Input file format not recognized, use .bam,.sam or .cram extension")

    click.echo("Input file: {}".format(args["bam"], kind), err=True)

    try:
        bam = pysam.AlignmentFile(args["bam"], opts[kind], threads=args["procs"])
    except:
        raise RuntimeError("Problem opening bam/sam/cram, check file has .bam/.sam/.cram in file name")

    bam_i = iter_bam(bam, args["search"])

    exc_tree = None
    if args["exclude"]:
        click.echo("Excluding {} from search".format(args["exclude"]), err=True)
        exc_tree = io_funcs.overlap_regions(args["exclude"])

    clip_length = args["clip_length"]

    send_output = None
    if args["output"] != "None":
        if args["output"] == "stdout":
            v = "-"
        else:
            v = args["output"]
        send_output = pysam.AlignmentFile(v, "wb", template=bam)

    out_name = "-" if (args["reads"] == "-" or args["reads"] == "stdout") else args["reads"]


    reads_out = pysam.AlignmentFile(out_name, "wbu", template=bam)

    # outbam = pysam.AlignmentFile(out_name + ".bam", "wbu", template=bam, threads=args["procs"])

    return bam, bam_i, exc_tree, clip_length, send_output, reads_out


cdef tuple get_reads(args):

    cdef uint32_t clip_length

    cdef int flag, tmpl, ct2, good_clip, nn
    cdef str qname
    cdef list cigartuples

    # Borrowed from lumpy
    cdef int required = 97
    restricted = 3484
    cdef int flag_mask = required | restricted

    bam, bam_i, exc_tree, clip_length, send_output, outbam = config(args)

    scope = deque([])
    read_names = set([])
    insert_size = []
    read_length = []
    cdef int max_scope = 100000

    cdef int count = 0
    bad = []
    for nn, r in enumerate(bam_i):

        if send_output:
            send_output.write(r)

        while len(scope) > max_scope:
            query = scope.popleft()
            if query.qname in read_names:
                outbam.write(query)
                count += 1

        scope.append(r)

        if exc_tree:  # Skip exclude regions
            if io_funcs.intersecter_int_chrom(exc_tree, r.rname, r.pos, r.pos + 1):
                continue

        flag = r.flag
        if flag & 1280:
            continue  # Read is duplicate or not primary alignment

        if len(insert_size) < 100000:

            if r.seq is not None:
                if r.rname == r.rnext and flag & flag_mask == required and r.tlen >= 0:
                    read_length.append(r.infer_read_length())
                    insert_size.append(r.tlen)

        qname = r.qname

        # if qname not in read_names:
        if not flag & 2 or flag & 2048:  # Save if read is discordant or supplementary
            read_names.add(qname)

        elif map_set_utils.cigar_clip(r, clip_length):
            read_names.add(qname)

        elif r.has_tag("SA"):
            read_names.add(qname)

    while len(scope) > 0:
        query = scope.popleft()
        if query.qname in read_names:
            outbam.write(query)
            count += 1

    outbam.close()
    if send_output:
        send_output.close()

    if len(insert_size) == 0:
        insert_median, insert_stdev = args["insert_median"], args["insert_stdev"]
        click.echo("WARNING: could not infer insert size, no 'normal' pairings found. Using arbitrary values", err=True)
    else:
        insert_median, insert_stdev = get_insert_params(insert_size)
        insert_median, insert_stdev = np.round(insert_median, 2), np.round(insert_stdev, 2)

    approx_read_length = int(np.mean(read_length))
    click.echo("Total input reads in bam {}".format(nn + 1), err=True)
    click.echo(f"Inferred read length {approx_read_length}, "
                   f"insert median {insert_median}, "
                   f"insert stdev {insert_stdev}", err=True)

    return count, insert_median, insert_stdev, approx_read_length


def process(args):
    t0 = time.time()

    count, insert_median, insert_stdev, read_length = get_reads(args)
    if args["index"] == "True" and args["reads"] not in "-,stdout":
        call("samtools index -@{} {}".format(args["procs"], args["reads"]), shell=True)

    click.echo("dysgu fetch {}, n={}, mem={} Mb, time={} h:m:s".format(args["bam"],
                                                                    count,
                                                                   int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
                                                                   str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)

    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            }

