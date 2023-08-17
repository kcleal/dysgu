#cython: language_level=3

from __future__ import absolute_import
import pysam
import time
import datetime
from collections import defaultdict
import os
import sys
import itertools
import logging
from libc.stdint cimport uint32_t

import pkg_resources
from dysgu.map_set_utils import echo
from dysgu.coverage import auto_max_cov
from dysgu.io_funcs import bed_iter


# thanks to senderle https://stackoverflow.com/questions/6462272/subtract-overlaps-between-two-ranges-without-sets
def range_diff(r1, r2):
    s1, e1 = r1
    s2, e2 = r2
    endpoints = sorted((s1, s2, e1, e2))
    result = []
    if endpoints[0] == s1 and endpoints[1] != s1:
        result.append((endpoints[0], endpoints[1]))
    if endpoints[3] == e1 and endpoints[2] != e1:
        result.append((endpoints[2], endpoints[3]))
    return result


def multirange_diff(r1_list, r2_list):
    for r2 in r2_list:
        r1_list = list(itertools.chain(*[range_diff(r1, r2) for r1 in r1_list]))
    return r1_list


def merge_simple(intervals):
    sorted_by_lower_bound = sorted(intervals)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def parse_search_regions(search, exclude, bam, first_delim=":", sep=","):

    is_chr = False
    chr_in_rname = any('chr' in i for i in bam.references)
    if search is not None and exclude is None:
        s = ""
        for line in bed_iter(search):
            if line[0] == "#":
                continue
            chrom, start, end = line.strip().split("\t", 4)[:3]
            if 'chr' in chrom:
                is_chr = True
            s += f"{chrom}{first_delim}{start}-{end}{sep}"
        if len(s) == 0:
            raise ValueError("Search regions not understood")
        return s

    targets = defaultdict(list)  # start end coords for chromosomes
    if search is not None:
        for line in bed_iter(search):
            if line[0] == "#":
                continue
            chrom, start, end = line.strip().split("\t", 4)[:3]
            if 'chr' in chrom:
                is_chr = True
            start = int(start)
            end = int(end)
            assert end > start
            targets[chrom].append((start, end))
    else:
        for l, c in zip(bam.lengths, bam.references):
            targets[c].append((0, l + 1))

    excl = defaultdict(list)
    if exclude is not None:
        for line in bed_iter(exclude):
        # with open(exclude, "r") as bed:
        #     for line in bed:
            if line[0] == "#":
                continue

            chrom, start, end = line.strip().split("\t", 4)[:3]
            if 'chr' in chrom:
                is_chr = True
            start = int(start)
            end = int(end)
            assert end > start
            excl[chrom].append((start, end))

    if chr_in_rname != is_chr:
        logging.warning('"chr" name conflict, please check --search/--exclude and bam chromosome names match')
    targets = {k: merge_simple(v) for k, v in targets.items()}
    excl = {k: merge_simple(v) for k, v in excl.items()}
    s = ""
    for chrom in targets:
        v = targets[chrom]
        if chrom in excl:
            v = multirange_diff(v, excl[chrom])

        for start, end in v:
            s += f"{chrom}{first_delim}{start}-{end}{sep}"

    if len(s) == 0:
        raise ValueError("Search/exclude regions not understood")
    return s


def assert_indexed_input(bam, fasta):
    if bam == "-" or bam == "stdin":
        raise ValueError("Cannot use a stream with current input options (--search/--exclude)")
    kind = bam.split(".")[-1]
    opts = {"bam": "rb", "cram": "rc", "sam": "r", "-": "r", "stdin": "r"}
    if kind not in opts:
        raise ValueError("Input file format not recognized, use .bam or .cram extension")
    try:
        bam = pysam.AlignmentFile(bam, opts[kind], reference_filename=fasta)
    except:
        raise RuntimeError("Problem opening input file, check file has .bam/.sam/.cram in file name, and file has a header")
    if not bam.has_index():
        raise ValueError("Input file must be indexed when options --search/--exclude are used")
    return bam


cdef extern from "find_reads.hpp":
    cdef int search_hts_alignments(char* infile, char* outfile, uint32_t min_within_size, int clip_length, int mapq_thresh,
                                   int threads, int paired_end, char* temp_f, int max_coverage, char* region,
                                   char* max_cov_ignore, char *fasta, bint write_all, char* out_write_mode_b)

def process(args):

    t0 = time.time()
    temp_dir = args["working_directory"]
    assert os.path.exists(temp_dir)

    if args["search"]:
        logging.info("Searching regions from {}".format(args["search"]))

    if args["exclude"]:
        logging.info("Excluding {} from search".format(args["exclude"]))

    if not args["output"]:
        bname = os.path.splitext(os.path.basename(args["bam"]))[0]
        if bname == "-":
            bname = os.path.basename(temp_dir)
        out_name = "{}/{}.{}.bam".format(temp_dir, bname, args["pfix"])

    else:
        bname = "-"
        out_name = args["output"]

    with open(os.path.join(temp_dir, "fetch_command.txt"), "w") as info:
        info.write(f"# dysgu v{pkg_resources.require('dysgu')[0].version}\n")
        info.write(" ".join(sys.argv))

    # update max cov automatically if applicable
    args["max_cov"] = auto_max_cov(args["max_cov"], args["bam"])
    pe = int(args["pl"] == "pe")

    cdef bytes infile_string_b = args["bam"].encode("ascii")
    cdef bytes fasta_b

    if "reference" in args:
        fasta_b = args["reference"].encode("ascii")
    else:
        fasta_b = "".encode("ascii")

    cdef bytes outfile_string_b = out_name.encode("ascii")
    cdef bytes out_write_mode_b = args["compression"].encode("ascii")
    cdef bytes temp_f = temp_dir.encode("ascii")
    cdef bint write_all = args["write_all"]
    cdef bytes regionbytes

    region = ".,"
    if args["search"] or args["exclude"]:
        bam = assert_indexed_input(args["bam"], args["reference"])
        region = parse_search_regions(args["search"], args["exclude"], bam)  # target regions in string format i.e. {chrom}:{start}-{end},

    regionbytes = region.encode("ascii")

    max_cov_ignore = ".,"
    max_cov_ignore_bytes = max_cov_ignore.encode("ascii")

    count = search_hts_alignments(infile_string_b,
                                  outfile_string_b,
                                  args["min_size"],
                                  args["clip_length"],
                                  args["mq"],
                                  args["procs"],
                                  pe,
                                  temp_f,
                                  int(args["max_cov"]),
                                  regionbytes,
                                  max_cov_ignore_bytes,
                                  fasta_b,
                                  write_all,
                                  out_write_mode_b)

    if count < 0:
        logging.critical("Error reading from input file, exit code {}".format(count))
        quit()
    elif count == 0:
        logging.critical("No reads found")
        quit()
    logging.info("dysgu fetch {} written to {}, n={}, time={} h:m:s".format(args["bam"], out_name,
                                                            count,
                                                            str(datetime.timedelta(seconds=int(time.time() - t0)))))
    return args["max_cov"]
