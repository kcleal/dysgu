#cython: language_level=3

from __future__ import absolute_import
from collections import defaultdict
import os
import click
import ncls
import numpy as np
from dysgu import samclips


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory {}".format(d))


cpdef dict make_template(rows, max_d, last_seen_chrom, fq, pairing_params, paired_end, isize, match_score, bias,
                         replace_hard):
    # Make a pickle-able data object for multiprocessing
    return {"isize": isize,
            "max_d": max_d,
            "match_score": match_score,
            "pairing_params": pairing_params,
            "paired_end": paired_end,
            "inputdata": rows,
            "bias": bias,
            "read1_length": 0,
            "read2_length": 0,
            "score_mat": {},
            "passed": 0,
            "name": rows[0][0][0],
            "last_seen_chrom": last_seen_chrom,
            "inputfq": fq,
            "read1_seq": 0,  # Some sam records may have seq == '*' , need a record of full seq for adding back in
            "read2_seq": 0,
            "read1_q": 0,
            "read2_q": 0,
            "read1_reverse": 0,  # Set to true if aligner has reverse complemented the sequence
            "read2_reverse": 0,
            "replace_hard": replace_hard,
            "fq_read1_seq": 0,
            "fq_read2_seq": 0,
            "fq_read1_q": 0,
            "fq_read2_q": 0,
            "read1_unmapped": 0,
            "read2_unmapped": 0
            }


cpdef list to_output(dict template):

    if "outstr" in template:
        return list(template["outstr"])

    return samclips.fixsam(template)


def sam_to_str(str template_name, list sam):
    return "".join(template_name + "\t" + "\t".join(i) + "\n" for i in sam)


def get_bed_regions(bed):
    b = [tuple([int(j) if j.isdigit() else j for j in i.strip().split("\t")[:3]]) for i in open(bed, "r")
         if i[0] != "#" and len(i) > 0 and "\t" in i]
    if len(b) == 0:
        raise ValueError("Bed regions not formatted correctly")
    return b


def overlap_regions(bed):
    if not bed:
        return None
    regions = get_bed_regions(bed)
    chrom_interval_start = defaultdict(list)
    chrom_interval_end = defaultdict(list)
    for c, s, e in regions:
        chrom_interval_start[c].append(int(s))
        chrom_interval_end[c].append(int(e))

    regions = {k: ncls.NCLS(np.array(chrom_interval_start[k]),
                            np.array(chrom_interval_end[k]),
                            np.array(chrom_interval_start[k])) for k in chrom_interval_start}

    return regions


def intersecter(tree, chrom, start, end):
    if tree is None:
        return 0
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return 1
        else:
            return 0
    else:
        return 0


def get_include_reads(include_regions, bam):

    if not include_regions:
        for r in bam:
            yield r

    regions = [i.strip().split("\t")[:3] for i in open(include_regions, "r") if i[0] != "#"]
    for c, s, e in regions:
        click.echo("Reading {}:{}-{}".format(c, s, e), err=True)
        for r in bam.fetch(c, int(s), int(e)):
            yield r


def sam_itr(args):

    itr = args["sam"]
    tree = overlap_regions(args["include"])

    # First get header
    header_string = ""
    last_seen_chrom = ""
    first_line = ""
    for t in itr:

        # t = str(t.decode("ascii"))

        if t[0] == "@":
            header_string += t
            continue

        first_line = t.split("\t", 4)
        last_seen_chrom = first_line[2]

        yield header_string
        break

    try:
        pos = int(first_line[3])
    except IndexError:
        raise IOError("No sam lines detected")

    ol = intersecter(tree, first_line[2], pos, pos + 250)

    yield first_line, last_seen_chrom, ol

    for t in itr:
        # t = str(t.decode("ascii"))
        line = t.split("\t", 4)

        if line[3] != last_seen_chrom:
            last_seen_chrom = line[2]

        pos = int(line[3])
        ol = intersecter(tree, line[2], pos, pos + 250)

        yield line, last_seen_chrom, ol


def fq_reader(args):
    # Iterate the fq files, send back a generator to use, generate only the read name, seq and qual lines
    # If > is at the start, iterate fasta instead
    def readfq(f):
        for l1 in f:
            l1 = l1.strip()
            last2 = l1[-2:]
            if l1[0] == "@":
                fastq = True
            else:
                fastq = False
            if last2 == "/2" or last2 == "/1":
                l1 = l1[1:-2]  # Strip trailing /1 or /2 and leading @
            else:
                l1 = l1[1:]
            l2 = next(f).strip()
            if fastq:
                next(f)  # Skip "+"
                l4 = next(f).strip()
            else:
                l4 = "1" * len(l2)
            yield l1, l2, l4

    if args["fq1"] is None:
        yield None, None

    if args["fq1"] and args["fq2"] is None:  # Single end
        with open(args["fq1"]) as fq1:
            for item in readfq(fq1):
                yield item, None
    else:
        with open(args["fq1"]) as fq1, open(args["fq2"]) as fq2:
            for item1, item2 in zip(readfq(fq1), readfq(fq2)):
                assert item1[0] == item2[0]
                yield item1, item2


def fq_getter(reader, name, args, fbuffer):

    if args["fq1"] is None:
        return None, None

    if name in fbuffer:
        fqlines = fbuffer[name]
        del fbuffer[fqlines[name]]
        return fqlines

    while True:
        fqlines = next(reader)
        q = fqlines[0][0]
        if q == name:
            return fqlines
        else:
            fbuffer[q] = fqlines


def iterate_mappings(args, version):

    params = {'clip_length', 'search', 'include', 'paired', 'max_insertion', 'min_aln', 'max_overlap',
    'ins_cost', 'ol_cost', 'inter_cost', 'u', 'match_score', 'bias', 'replace_hardclips', 'fq1', 'fq2',
     'insert_median', 'insert_stdev', 'mq', 'max_tlen', 'template_size'}
    cp_args = {k: v for k, v in args.items() if k in params}

    arg_str = ", ".join(["{}={}".format(i, j) for i, j in args.items() if i in params])
    inputstream = sam_itr(args)

    total = 0
    name = ""
    rows = []
    header_string = next(inputstream)
    header_string += "@PG\tID:DYSGU\tPN:dysgu choose\tVN:{}\tCL:{}\n".format(version, arg_str)

    yield header_string

    fq_buffer = defaultdict(list)
    fq_iter = fq_reader(args)

    last_seen_chrom = ""

    for m, last_seen_chrom, ol in inputstream:  # Alignment

        nm = m[0]
        if name != nm:
            if len(rows) > 0:
                total += 1
                fq = fq_getter(fq_iter, name, args, fq_buffer)

                yield rows, last_seen_chrom, fq

            rows = []
            name = nm

        rows.append((m, ol))  # String, ol states if alignment overlaps ROI

    # Deal with last record
    if len(rows) > 0:
        total += 1
        fq = fq_getter(fq_iter, name, args, fq_buffer)
        yield rows, last_seen_chrom, fq

    click.echo("Total processed " + str(total), err=True)
