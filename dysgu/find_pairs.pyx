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

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment

from libc.stdint cimport uint32_t

from dysgu import data_io, io_funcs
from dysgu cimport map_set_utils


# def iter_bam(AlignmentFile bam, search, show=True):
def iter_bam(bam, search, show=True):
    # cdef AlignedSegment aln

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
    click.echo("Input file: {}".format(args["bam"], kind), err=True)

    # cdef AlignmentFile bam

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
    if args["dest"]:

        if "/" in args["bam"]:
            out_name = args["dest"] + "/" + args["bam"].rsplit("/", 1)[1].rsplit(".", 1)[0] + "." + args["post_fix"]
        else:
            out_name = args["dest"] + "/" + args["bam"].rsplit(".", 1)[0] + "." + args["post_fix"]
        if not os.path.exists(args["dest"]):
            os.makedirs(args["dest"])
    else:
        if "/" in args["bam"]:
            out_name = os.getcwd() + "/" + args["bam"].rsplit("/", 1)[1].rsplit(".", 1)[0] + "." + args["post_fix"]
        else:
            out_name = os.getcwd() + "/" + args["bam"].rsplit(".", 1)[0] + "." + args["post_fix"]

    if args["mapper"] == "last":
        call("samtools view -H -o {}.dict {}".format(out_name, args["bam"]), shell=True)

    outbam = pysam.AlignmentFile(out_name + ".bam", "wbu", template=bam, threads=args["procs"])

    return bam, bam_i, exc_tree, clip_length, out_name, outbam


cdef tuple get_reads(args):

    cdef uint32_t clip_length
    cdef str out_name
    # cdef AlignmentFile bam, outbam
    cdef int flag, tmpl, ct2, good_clip
    cdef str qname
    cdef list cigartuples
    # cdef AlignedSegment r

    bam, bam_i, exc_tree, clip_length, out_name, outbam = config(args)

    scope = deque([])
    read_names = set([])
    insert_size = []
    read_length = []
    cdef int max_scope = 20000

    for r in bam_i:

        if exc_tree:  # Skip exclude regions
            if io_funcs.intersecter_int_chrom(exc_tree, r.rname, r.pos, r.pos + 1):
                continue

        flag = r.flag
        if flag & 1280:
            continue  # Read is duplicate or not primary alignment

        if len(insert_size) < 100000:
            if flag & 2:
                tmpl = abs(r.template_length)
                if tmpl < 2000:
                    insert_size.append(tmpl)
        if len(read_length) < 100000:
            rl = r.infer_read_length()
            if rl:
                read_length.append(rl)

        qname = r.qname

        # if qname not in read_names:
        if not flag & 2 or flag & 2048:  # Save if read is discordant or supplementary
                read_names.add(qname)

        elif map_set_utils.cigar_clip(r, clip_length):
                read_names.add(qname)

        elif r.has_tag("SA"):
            read_names.add(qname)

        scope.append(r)

        while len(scope) > max_scope:
            query = scope.popleft()
            if query.qname in read_names:
                outbam.write(query)

    while len(scope) > 0:
        query = scope.popleft()
        if query.qname in read_names:
            outbam.write(query)

    outbam.close()

    if len(insert_size) == 0:
        insert_median, insert_stdev = args["insert_median"], args["insert_stdev"]
        click.echo("WARNING: could not infer insert size, no 'normal' pairings found. Using arbitrary values", err=True)
    else:
        insert_median, insert_stdev = np.median(insert_size), np.std(insert_size)
        insert_median, insert_stdev = np.round(insert_median, 2), np.round(insert_stdev, 2)

    click.echo("Median insert size: {} (+/- {})".format(insert_median, insert_stdev), err=True)

    return insert_median, insert_stdev, int(np.mean(read_length)), out_name


def process(args):
    t0 = time.time()

    data_io.mk_dest(args["dest"])

    insert_median, insert_stdev, read_length, out_name = get_reads(args)

    call("samtools index -@{} {}.bam".format(args["procs"], out_name), shell=True)

    click.echo("Collected reads in {} h:m:s".format(str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)

    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            "fastq": out_name,
            "out_pfix": out_name}

