#cython: language_level = 3
"""
Compiled utils to generate proper sam output and flag information
"""
from __future__ import absolute_import

import pysam
import time
import datetime
import numpy as np
import os
import click
from subprocess import call
from dysgu import data_io, io_funcs

# from pysam.libcalignedsegment cimport AlignedSegment


def echo(*args):
    click.echo(args, err=True)


def iter_bam(bam, search):
    if not search:
        click.echo("Searching input file", err=True)
        for aln in bam.fetch(until_eof=True):  # Also get unmapped reads
            yield aln
    else:
        click.echo("Limiting search to {}".format(search), err=True)
        with open(search, "r") as bed:
            for line in bed:
                if line[0] == "#":
                    continue
                chrom, start, end = line.split("\t")[:3]
                for aln in bam.fetch(reference=chrom, start=int(start), end=int(end)):  # Also get unmapped reads
                    yield aln


def get_reads_f(args):

    paired = args["paired"] == "True"
    fq_out = args["out_format"] == "fq"

    kind = args["bam"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "r"}
    click.echo("Input file is {}, (.{} format)".format(args["bam"], kind), err=True)

    bam = pysam.AlignmentFile(args["bam"], opts[kind])
    bam_i = iter_bam(bam, args["search"])

    if args["exclude"]:
        click.echo("Excluding {} from search".format(args["exclude"]), err=True)
        exc_tree = data_io.overlap_regions(args["exclude"])


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


    def container(paired):
        if paired:
            # write_me, alignments
            return [0, None, None]
        return [0, None]

    cdef dict read_cache = {}

    cdef list insert_size = []
    cdef list read_length = []

    cdef int bam_only = 0
    cdef int clip_length = args["clip_length"]
    cdef int tmpl, flag, pos

    cdef int aligns_written = 0
    cdef str qname = ""

    cdef list d, ct

    if paired:
        fq1 = open(out_name + f"1.{args['out_format']}", "w")
        fq2 = open(out_name + f"2.{args['out_format']}", "w")
    else:
        fq1 = open(out_name + f".{args['out_format']}", "w")

    for r in bam_i:

        if args["exclude"]:  # Skip exclude regions
            pos = r.pos
            if data_io.intersecter(exc_tree, bam.get_reference_name(r.rname), pos, pos + 1):
                continue

        flag = r.flag

        if flag & 3328 or not r.cigar: #r.cigartuples:
            continue  # Read is duplicate or not primary alignment or supplementary, Cigar is not formatted

        if paired and len(insert_size) < 10000:
            if flag & 2:
                tmpl = abs(r.template_length)
                if tmpl < 1000:
                    insert_size.append(tmpl)

            read_length.append(r.infer_read_length())

        qname = r.qname + str(pos) + str(flag)

        if qname not in read_cache:
            d = container(paired)

        else:
            d = read_cache[qname]

        if flag & 64 or not paired:
            d[1] = r

        elif paired:
            d[2] = r

        # Check for discordants
        if not flag & 2:
            d[0] = 1

        else:
            # Check for soft-clips
            ct = r.cigartuples
            if (ct[0][0] == 4 and ct[0][1] >= clip_length) or (ct[-1][0] == 4 and ct[-1][1] >= clip_length):
                d[0] = 1

            # Check for SA tag
            elif r.has_tag("SA"):
                d[0] = 1

        if (paired and d[1] is not None and d[2] is not None) or (not paired and d[1] is not None):
            if d[0]:
                s1 = d[1].seq
                q1 = d[1].qual
                q_exists = True
                if q1 is None and fq_out:
                    q_exists = False
                    q1 = "1" * len(s1)
                if d[1].flag & 16:
                    s1 = io_funcs.reverse_complement(s1, len(s1))
                    if q_exists and fq_out:
                        q1 = q1[::-1]

                if fq_out:
                    fq1.write(f"@{qname}\n{s1}\n+\n{q1}\n")
                else:
                    fq1.write(f">{qname}\n{s1}\n")

                if not paired:
                    aligns_written += 1
                    continue

                s2 = d[2].seq
                q2 = d[2].qual
                q2_exists = True
                if q2 is None and fq_out:
                    q2_exists = False
                    q2 = "1" * len(s1)
                if d[2].flag & 16:
                    s2 = io_funcs.reverse_complement(s2, len(s2))
                    if q2_exists and fq_out:
                        q2 = q2[::-1]

                if fq_out:
                    fq2.write(f"@{qname}\n{s2}\n+\n{q2}\n")
                else:
                    fq2.write(f">{qname}\n{s2}\n")

                aligns_written += 2
            if paired:
                del read_cache[qname]

        elif paired:
            read_cache[qname] = d

    if aligns_written == 0:
        click.echo("No reads found, aborting", err=True)
        quit()
    click.echo("Writen {} alignments".format(aligns_written), err=True)

    fq1.close()
    if paired:
        fq2.close()
    bam.close()

    if paired:
        if len(insert_size) == 0:
            insert_median, insert_stdev = args["insert_median"], args["insert_stdev"]
            click.echo("WARNING: could not infer insert size, no 'normal' pairings found. Using arbitrary values", err=True)
        else:
            insert_median, insert_stdev = np.median(insert_size), np.std(insert_size)
            insert_median, insert_stdev = np.round(insert_median, 2), np.round(insert_stdev, 2)

        click.echo("Median insert size: {} (+/- {})".format(insert_median, insert_stdev), err=True)

        return insert_median, insert_stdev, int(np.mean(read_length)), out_name

    else:
        return 0, 0, 0, out_name


def process(args):

    t0 = time.time()
    data_io.mk_dest(args["dest"])
    insert_median, insert_stdev, read_length, out_name = get_reads_f(args)
    click.echo("Collected reads in {} h:m:s".format(str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)

    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            "fastq": out_name,
            "out_pfix": out_name}
