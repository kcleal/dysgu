#cython: language_level = 3
"""
Compiled utils to generate proper sam output and flag information
"""
from __future__ import absolute_import
import pysam
import time
import datetime
import numpy as np
import click
import resource
from subprocess import run
from sys import stdout
from dysgu import io_funcs
from dysgu.coverage import get_insert_params


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
    # todo bug sometimes duplicate fq reads occur?
    two_files = True if args["reads2"] != "None" else False

    fq_out = args["out_format"] == "fq"

    kind = args["bam"].split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "r"}
    if kind not in opts:
        raise ValueError("Input file format not recognized, use .bam,.sam or .cram extension")
    click.echo("Input file: {}, (.{} format). Reading threads={}".format(args["bam"], kind, args["procs"]), err=True)

    bam = pysam.AlignmentFile(args["bam"], opts[kind], threads=args["procs"])
    bam_i = iter_bam(bam, args["search"])

    exc_tree = None
    if args["exclude"]:
        click.echo("Excluding {} from search".format(args["exclude"]), err=True)
        exc_tree = io_funcs.overlap_regionspy(args["exclude"])


    def container(paired_end):
        if paired_end:
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
    cdef int nn = 0

    if two_files:
        assert f"{args['reads']}" != f"{args['reads2']}"
        fq1 = open(f"{args['reads']}", "w")
        fq2 = open(f"{args['reads2']}", "w")
    else:
        if args["reads"] in "-stdout":
            fq1 = stdout
        else:
            fq1 = open(f"{args['reads']}", "w")  # Write interleaved or single end

    # Borrowed from lumpy
    cdef int required = 97
    restricted = 3484
    cdef int flag_mask = required | restricted

    for nn, r in enumerate(bam_i):

        if args["exclude"]:  # Skip exclude regions
            pos = r.pos
            if io_funcs.intersecterpy(exc_tree, bam.get_reference_name(r.rname), pos, pos + 1):
                continue

        flag = r.flag

        paired_end = False
        if flag & 131:  # read_paired, proper_pair, second in pair
            paired_end = True

        if flag & 3328 or not r.cigar:
            continue  # Read is duplicate or not primary alignment or supplementary, Cigar is not formatted

        if paired_end and len(insert_size) < 10000:
            if r.seq is not None:
                if r.rname == r.rnext and flag & flag_mask == required and r.tlen >= 0:
                    read_length.append(r.infer_read_length())
                    insert_size.append(r.tlen)

        qname = r.qname

        if qname not in read_cache:
            d = container(paired_end)

        else:
            d = read_cache[qname]

        if flag & 64 or not paired_end:
            d[1] = r

        elif paired_end:
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

        if (paired_end and d[1] is not None and d[2] is not None) or (not paired_end and d[1] is not None):
            if d[0]:  # --> write to output

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
                    fq1.write(f"@{qname}/1\n{s1}\n+\n{q1}\n")
                else:
                    fq1.write(f">{qname}/1\n{s1}\n")

                if not paired_end:
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
                    if two_files:
                        fq2.write(f"@{qname}/2\n{s2}\n+\n{q2}\n")
                    else:
                        fq1.write(f"@{qname}/2\n{s2}\n+\n{q2}\n")
                else:
                    if two_files:
                        fq2.write(f">{qname}/2\n{s2}\n")
                    else:
                        fq1.write(f">{qname}/2\n{s2}\n")

                aligns_written += 2
            #if paired_end:
            if qname in read_cache:  #  drop bad reads or written
                del read_cache[qname]  # Only cached if paired end

        elif paired_end:
            read_cache[qname] = d

    #

    if aligns_written == 0:
        click.echo("No reads found, aborting", err=True)
        quit()

    fq1.close()
    if two_files:
        fq2.close()
    bam.close()

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

    return insert_median, insert_stdev, approx_read_length, aligns_written


def process(args):

    t0 = time.time()
    io_funcs.mk_dest(args["dest"])
    click.echo("Output format: {}".format(args["out_format"]), err=True)
    insert_median, insert_stdev, read_length, count = get_reads_f(args)

    click.echo("dysgu fetch {}, n={}, mem={} Mb, time={} h:m:s".format(args["bam"],
                                                                    count,
                                                                   int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6),
                                                                   str(datetime.timedelta(seconds=int(time.time() - t0)))), err=True)

    return {"insert_median": insert_median,
            "insert_stdev": insert_stdev,
            "read_length": read_length,
            }
