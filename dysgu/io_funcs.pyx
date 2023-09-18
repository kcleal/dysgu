#cython: language_level=2, boundscheck=True, wraparound=True
#distutils: language=c++

import numpy as np
cimport numpy as np
import logging
from map_set_utils import merge_intervals, echo
from collections import defaultdict
from importlib.metadata import version
import sortedcontainers
import pandas as pd
import os
import sys
import gzip
from dysgu.map_set_utils import Py_BasicIntervalTree
import random

from libc.stdlib cimport malloc


cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0',
                       '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0',  'a',  'a' ]

np.random.seed(0)
random.seed(0)

cpdef str reverse_complement(str seq, int seq_len):
    """https://bioinformatics.stackexchange.com/questions/3583/\
    what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho/3595#3595"""

    cdef char *seq_dest = <char *>malloc(seq_len + 1)
    seq_dest[seq_len] = '\0'

    cdef bytes py_bytes = seq.encode('UTF-8')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('UTF-8')


def bed_iter(path):
    if path[-4:] == ".bed":
        with open(path, "r") as bed:
            for line in bed:
                yield line
    elif path[-3:] == ".gz":
        with gzip.open(path, "r") as bed:
            for line in bed:
                yield line.decode('utf-8')
    else:
        raise ValueError("File extension not understood. Input bed file must have extension .gz or .bed")


cpdef list get_bed_regions(bed):
    b = []
    for line in bed_iter(bed):
        if line[0] != "#":
            r = line.strip().split("\t")
            b.append((r[0], int(r[1]), int(r[2])))
    if len(b) == 0:
        raise ValueError("Bed regions empty")
    return b


def iitree(a, add_value=False):
    # sorted input list a will not lead to a balanced binary-search-tree if added in sequential order
    # This function reorders the list to ensure a balanced BST when added using Py_BasicIntervalTree
    tree = Py_BasicIntervalTree()
    index = 0
    for h in a:
        if add_value:
            tree.add(h[0], h[1], h[2])
        else:
            tree.add(h[0], h[1], index)
        index += 1
    tree.index()
    return tree


cpdef dict overlap_regions(str bed, int_chroms=False, infile=None):
    if not bed:
        return {}
    regions = merge_intervals(get_bed_regions(bed))
    chrom_intervals = defaultdict(list)
    for c, s, e in regions:
        if int_chroms:
            c = infile.gettid(c)
        chrom_intervals[c].append((s, e))

    # create balanced ordering
    chrom_intervals = {k: iitree(v) for k, v in chrom_intervals.items()}
    return chrom_intervals


cpdef int intersecter(tree, chrom, int start, int end):
    cdef bint found = 0
    if tree and chrom in tree:
        found = tree[chrom].searchInterval(start, end)
    return found


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory {}".format(d))


def get_include_reads(include_regions, bam):
    if not include_regions:
        for r in bam:
            yield r
    regions = [i.strip().split("\t")[:3] for i in open(include_regions, "r") if i[0] != "#"]
    for c, s, e in regions:
        logging.info("Reading {}:{}-{}".format(c, s, e))
        for r in bam.fetch(c, int(s), int(e)):
            yield r


cpdef list col_names(small_output):

    if small_output:
        return ["chrA", "posA", "chrB", "posB", "sample", "event_id", "grp_id", "n_in_grp", "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B", 'contigA', 'contigB', "svlen", "svlen_precise", "rep", "gc",
          ["GT", "GQ", "MAPQpri", "su", "spanning", "pe", "supp", "sc", "bnd",
         "raw_reads_10kb", "neigh10kb", "plus",
                "minus", "remap_score", "remap_ed", "bad_clip_count", "fcc", "inner_cn", "outer_cn", "prob"]
            ]
    else:
        return ["chrA", "posA", "chrB", "posB", "sample", "event_id", "grp_id", "n_in_grp", "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B", 'contigA', 'contigB', "svlen", "svlen_precise",  "rep", "gc",
          ["GT", "GQ", "NMpri", "NMsupp", "NMbase", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "su", "spanning", "pe", "supp", "sc", "bnd", "sqc", "scw", "clip_qual_ratio", "block_edge",
         "raw_reads_10kb", "mcov",
          "linked", "neigh", "neigh10kb",  "ref_bases", "plus",
                "minus", "n_gaps", "n_sa", "n_xa", "n_unmapped_mates", "double_clips", "remap_score", "remap_ed", "bad_clip_count", "fcc", "n_small_tlen", "ras", "fas",
                "inner_cn", "outer_cn", "compress", "ref_rep", "prob"]
            ]


def make_main_record(r, dysgu_version, index, format_f, df_rows, add_kind, small_output):
    rep, repsc, lenprec = 0, 0, 1
    mean_prob, max_prob = None, None
    debug = False
    # if abs(r['posA'] - 39084726) < 10:
    #     echo('found', dict(r))
    #     echo(len(format_f) > 1)
    #     echo([(int(v["su"]), v["posA"], v["event_id"], k) for k, v in df_rows.items()])
    #     debug = True
    if len(format_f) > 1:
        best = sorted([(int(v["su"]), k) for k, v in df_rows.items()], reverse=True)[0][1]
        probs = [v["prob"] for k, v in df_rows.items()]
        mean_prob = np.mean(probs)
        max_prob = np.max(probs)
        r = df_rows[best]
        gc = round(r["gc"], 2)
        if gc == -1:
            gc = "."
        if not small_output:
            rep = r["rep"]
            repsc = r["rep_sc"]
        lenprec = r["svlen_precise"]
        n_expansion = r["n_expansion"]
        stride = r["stride"]
        exp_seq = r["exp_seq"]
        ref_poly = r["ref_poly_bases"]
        overlaps = r["query_overlap"]
        su, pe, sr, sc, bnd, wr = 0, 0, 0, 0, 0, 0
        for row in df_rows.values():
            pe += row["pe"]
            sr += row["supp"]
            sc += row["sc"]
            su += row["su"]
            bnd += row["bnd"]
            wr += row["spanning"]

    else:
        pe = r["pe"]
        sr = r["supp"]
        sc = r["sc"]
        su = r["su"]
        bnd = r["bnd"]
        wr = r["spanning"]
        gc = round(r["gc"], 2)
        if gc == -1:
            gc = "."
        if not small_output:
            rep = r["rep"]
            repsc = r["rep_sc"]
        lenprec = r["svlen_precise"]
        n_expansion = r["n_expansion"]
        stride = r["stride"]
        exp_seq = r["exp_seq"]
        ref_poly = r["ref_poly_bases"]
        overlaps = r["query_overlap"]

    samp = r["sample"]
    read_kind = r["type"]

    r["posA"] = max(1, r["posA"] + 1)  # convert to 1 based indexing
    r["posB"] = max(1, r["posB"] + 1)

    info_extras = []
    if r["chrA"] == r["chrB"]:
        if 'svlen' in r:
            info_extras.append(f"SVLEN={r['svlen']}")
        else:
            info_extras.append(f"SVLEN=0")

    if r["contigA"]:
        info_extras.append(f"CONTIGA={r['contigA']}")
    if r["contigB"]:
        info_extras.append(f"CONTIGB={r['contigB']}")

    if not r["variant_seq"] or r["variant_seq"][0] == "<":
        if "left_ins_seq" in r and r["left_ins_seq"]:
            info_extras.append(f"LEFT_SVINSSEQ={r['left_ins_seq']}")
        if "right_ins_seq" in r and r["right_ins_seq"]:
            info_extras.append(f"RIGHT_SVINSSEQ={r['right_ins_seq']}")

    if add_kind:
        info_extras += [f"KIND={r['kind']}"]

    if not small_output:

        try:
            info_extras.append(f"REP={'%.3f' % float(rep)}")
        except ValueError:  # rep is "."
            info_extras.append("REP=0")

        try:
            info_extras.append(f"REPSC={'%.3f' % float(repsc)}")
        except ValueError:
            info_extras.append("REPSC=0")

    info_extras += [
                    f"GC={gc}",
                    f"NEXP={n_expansion}",
                    f"STRIDE={stride}",
                    f"EXPSEQ={exp_seq}",
                    f"RPOLY={ref_poly}",
                    f"OL={overlaps}",
                    f"SU={su}",
                    f"WR={wr}",
                    f"PE={pe}",
                    f"SR={sr}",
                    f"SC={sc}",
                    f"BND={bnd}",
                    f"LPREC={lenprec}",
                    f"RT={read_kind}"]

    if r['chrA'] != r['chrB']:
        chr2_pos = f";CHR2_POS={r['posB']}"
    else:
        chr2_pos = ""

    if mean_prob is not None:
        info_extras += [f"MeanPROB={round(mean_prob, 3)}", f"MaxPROB={round(max_prob, 3)}"]

    if small_output:
        fmt_keys = "GT:GQ:MAPQP:SU:WR:PE:SR:SC:BND:COV:NEIGH10:PS:MS:RMS:RED:BCC:FCC:ICN:OCN:PROB"
    else:
        fmt_keys = "GT:GQ:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BND:SQC:SCW:SQR:BE:COV:MCOV:LNK:NEIGH:NEIGH10:RB:PS:MS:SBT:NG:NSA:NXA:NMU:NDC:RMS:RED:BCC:FCC:STL:RAS:FAS:ICN:OCN:CMP:RR:JIT:PROB"

    if "variant_seq" in r and isinstance(r["variant_seq"], str):
        if r['svtype'] == "INS" or r.variant_seq:
            alt_field = r.variant_seq.upper()
        else:
            alt_field = f"<{r['svtype']}>"
    else:
        alt_field = f"<{r['svtype']}>"

    if "ref_seq" in r and isinstance(r["ref_seq"], str):
        ref_field = r["ref_seq"]
    else:
        ref_field = "."

    rec = [r["chrA"], r["posA"], r["event_id"],
           ref_field,
           alt_field,
           ".", "." if "filter" not in r else r['filter'],
           # INFO line
           ";".join([f"SVMETHOD=DYSGUv{dysgu_version}",
                   f"SVTYPE={r['svtype']}",
                   f"END={r['posB']}" if r['chrA'] == r['chrB'] else f"END={r['posA'] + 1}",
                   f"CHR2={r['chrB']}" + chr2_pos,
                   f"GRP={r['grp_id']}",
                   f"NGRP={r['n_in_grp']}",
                   f"CT={r['join_type']}",
                   f"CIPOS95={r['cipos95A']}",
                   f"CIEND95={r['cipos95B']}",

                   ] + info_extras),
           fmt_keys  # :PROB
           ]
    # FORMAT line(s)
    for item in format_f.values():
        rec.append(":".join(map(str, item)))
    # if debug:
    #     echo("returned ", rec)
    return rec


def get_fmt(r, small_output):
    if small_output:
        v = [r["GT"], r["GQ"], r['MAPQpri'], r['su'], r['spanning'], r['pe'], r['supp'], r['sc'], r['bnd'],
             r['raw_reads_10kb'], r['neigh10kb'], r["plus"], r["minus"], r["remap_score"], r["remap_ed"],
             r["bad_clip_count"], round(r["fcc"], 3), round(r["inner_cn"], 3), round(r["outer_cn"], 3), r['prob']
             ]
        return v
    else:
        v = [r["GT"], r["GQ"], r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
             r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
             r['sc'], r['bnd'], round(r['sqc'], 2), round(r['scw'], 1), round(r['clip_qual_ratio'], 3), r['block_edge'], r['raw_reads_10kb'], round(r['mcov'], 2), int(r['linked']), r['neigh'], r['neigh10kb'],
             r['ref_bases'], r["plus"], r["minus"], round(r["strand_binom_t"], 4), r['n_gaps'], round(r["n_sa"], 2),
             round(r["n_xa"], 2), round(r["n_unmapped_mates"], 2), r["double_clips"], r["remap_score"], r["remap_ed"], r["bad_clip_count"], round(r["fcc"], 3), r["n_small_tlen"], r["ras"], r['fas'],
             round(r["inner_cn"], 3), round(r["outer_cn"], 3), round(r["compress"], 2), round(r["ref_rep"], 3), round(r["jitter"], 3), r['prob']
             ]
        return v


def gen_format_fields(r, df, names, n_fields, small_output):
    if len(names) == 1:
        return {0: get_fmt(r, small_output)}, {}
    cols = {}
    if "partners" in r:
        if not isinstance(r["partners"], set):
            if len(r["partners"]) == 0 or pd.isna(r["partners"]):
                r["partners"] = []
            else:
                r["partners"] = [int(i.split(",")[1]) for i in r["partners"].split("|")]
        for idx in r["partners"]:
            if idx in df.index:  # might be already dropped
                r2 = df.loc[idx]
                cols[r2["table_name"]] = r2
    if "table_name" in r:
        cols[r["table_name"]] = r
    format_fields = sortedcontainers.SortedDict()
    for name in names:
        if name in cols:
            row = cols[name]
            format_fields[name] = get_fmt(row, small_output)
        else:
            format_fields[name] = ['0/0'] + [0] * (n_fields - 1)
    return format_fields, cols


def get_header():
    header = """##fileformat=VCFv4.2
##source=DYSGU
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=lowProb,Description="Probability below threshold set with --thresholds">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=CHR2_POS,Number=1,Type=Integer,Description="Chromosome position for END coordinate in case of a translocation">
##INFO=<ID=GRP,Number=1,Type=Integer,Description="Group id for complex SVs">
##INFO=<ID=NGRP,Number=1,Type=Integer,Description="Number of SVs in group">
##INFO=<ID=RT,Number=1,Type=String,Description="Type of input reads, 1=paired-end, 2=pacbio, 3=nanopore">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=CIPOS95,Number=1,Type=Integer,Description="Confidence interval size (95%) around POS for imprecise variants">
##INFO=<ID=CIEND95,Number=1,Type=Integer,Description="Confidence interval size (95%) around END for imprecise variants">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=KIND,Number=1,Type=String,Description="Kind of join with respect to input regions">
##INFO=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=WR,Number=1,Type=Integer,Description="Number of reads that have SV within-read">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of supplementary reads supporting the variant across all samples">
##INFO=<ID=SC,Number=1,Type=Integer,Description="Number of soft-clip reads supporting the variant across all samples">
##INFO=<ID=BND,Number=1,Type=Integer,Description="Number of break-end alignments supporting the variant">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Contig from CHROM POS">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Contig from CHR2 END">
##INFO=<ID=GC,Number=1,Type=Float,Description="GC% of assembled contigs">
##INFO=<ID=REP,Number=1,Type=Float,Description="Repeat score for contigs aligned bases">
##INFO=<ID=REPSC,Number=1,Type=Float,Description="Repeat score for contigs soft-clipped bases">
##INFO=<ID=LPREC,Number=1,Type=Integer,Description="SV length precise=1, inferred=0">
##INFO=<ID=NEXP,Number=1,Type=Integer,Description="Number of expanded repeat bases at break">
##INFO=<ID=STRIDE,Number=1,Type=Integer,Description="Repeat expansion stride or period">
##INFO=<ID=EXPSEQ,Number=1,Type=String,Description="Expansion sequence">
##INFO=<ID=RPOLY,Number=1,Type=Integer,Description="Number of reference polymer bases">
##INFO=<ID=OL,Number=1,Type=Integer,Description="Query overlap in bp">
##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion">
##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description="Known left side of insertion for an insertion of unknown length">
##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description="Known right side of insertion for an insertion of unknown length">
##INFO=<ID=MeanPROB,Number=1,Type=Float,Description="Mean probability of event being true across samples">
##INFO=<ID=MaxPROB,Number=1,Type=Float,Description="Max probability of event being true across samples">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality phred scaled">
##FORMAT=<ID=NMP,Number=1,Type=Float,Description="Mean edit distance for primary alignments supporting the variant">
##FORMAT=<ID=NMS,Number=1,Type=Float,Description="Mean edit distance for supplementary alignments supporting the variant">
##FORMAT=<ID=NMB,Number=1,Type=Float,Description="Mean basic, edit distance. Gaps >= 30 bp are ignored">
##FORMAT=<ID=MAPQP,Number=1,Type=Float,Description="Mean MAPQ for primary reads supporting the variant">
##FORMAT=<ID=MAPQS,Number=1,Type=Float,Description="Mean MAPQ for supplementary reads supporting the variant">
##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Number of alignments in normal-pair orientation supporting the variant">
##FORMAT=<ID=MAS,Number=1,Type=Integer,Description="Maximum alignment score of supplementary reads supporting the variant">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=WR,Number=1,Type=Integer,Description="Number of reads that have SV within-read">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supplementary alignments supporting the variant">
##FORMAT=<ID=SC,Number=1,Type=Integer,Description="Number of soft-clipped alignments supporting the variant">
##FORMAT=<ID=BND,Number=1,Type=Integer,Description="Number of break-end alignments supporting the variant">
##FORMAT=<ID=SQC,Number=1,Type=Float,Description="Soft-clip quality value correlation between reads">
##FORMAT=<ID=SCW,Number=1,Type=Float,Description="Soft-clip quality weight value">
##FORMAT=<ID=SQR,Number=1,Type=Float,Description="Soft-clip base-quality ratio wrt to aligned bases">
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Mean read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=MCOV,Number=1,Type=Float,Description="Maximum read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other break points within 1 bp of break site">
##FORMAT=<ID=NEIGH10,Number=1,Type=Integer,Description="Number of other break points within 10 kp of break site">
##FORMAT=<ID=RB,Number=1,Type=Integer,Description="Number of reference bases in contigs">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Number of reads on plus strand">
##FORMAT=<ID=MS,Number=1,Type=Integer,Description="Number of reads on minus strand">
##FORMAT=<ID=SBT,Number=1,Type=Float,Description="Strand-bias, one-sided binomial test p-value">
##FORMAT=<ID=NG,Number=1,Type=Float,Description="Mean number of small gaps < 30 bp">
##FORMAT=<ID=NSA,Number=1,Type=Float,Description="Mean number of SA tags per read">
##FORMAT=<ID=NXA,Number=1,Type=Float,Description="Mean number of XA tags per read">
##FORMAT=<ID=NMU,Number=1,Type=Float,Description="Mean number of mates unmapped per read">
##FORMAT=<ID=NDC,Number=1,Type=Integer,Description="Number of double-clips, alignments with left and right clips">
##FORMAT=<ID=RMS,Number=1,Type=Integer,Description="Remapping score">
##FORMAT=<ID=RED,Number=1,Type=Integer,Description="Remapping edit distance">
##FORMAT=<ID=BCC,Number=1,Type=Integer,Description="Bad soft-clip count within +/- 500 bp">
##FORMAT=<ID=FCC,Number=1,Type=Float,Description="Fold-coverage change for SVs">
##FORMAT=<ID=STL,Number=1,Type=Integer,Description="N reads with small TLEN below 0.05% of distribution">
##FORMAT=<ID=RAS,Number=1,Type=Integer,Description="Reverse soft-clip to alignment score">
##FORMAT=<ID=FAS,Number=1,Type=Integer,Description="Forward soft-clip to alignment score">
##FORMAT=<ID=ICN,Number=1,Type=Float,Description="Inner copy number">
##FORMAT=<ID=OCN,Number=1,Type=Float,Description="Outer copy number">
##FORMAT=<ID=CMP,Number=1,Type=Float,Description="Compression ratio of contigs">
##FORMAT=<ID=RR,Number=1,Type=Float,Description="Repeat score for reference">
##FORMAT=<ID=JIT,Number=1,Type=Float,Description="SV length jitter">
##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event being true">{contig_names}
##command="{input_command}"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

    return header


def to_vcf(df, args, names, outfile, show_names=True,  contig_names="", header=None,
           small_output_f=True, sort_output=True, n_fields=None):

    header = get_header()
    input_command = ' '.join(sys.argv)
    outfile.write(
        header.format(input_command=input_command,
                      contig_names=contig_names) + "\t" + "\t".join(names) + "\n")

    if show_names:
        logging.info("Samples: {}".format(str(list(names))))
    dysgu_version = version("dysgu")
    seen_idx = set([])
    cnames = ['raw_reads_10kb', 'NMpri', 'NMsupp', 'MAPQpri', 'MAPQsupp', "NMbase", "n_gaps"]
    for col in cnames:
        if col in df.columns:
            df[col] = df[col].astype(float).round(2)

    for col in ['maxASsupp', 'neigh', 'neigh10kb']:
        if col in df.columns:
            df[col] = [int(i) if (i == i and i is not None) else 0 for i in df[col]]

    count = 0
    recs = []
    if args is not None:
        add_kind = args["add_kind"] == "True"
        if args["metrics"]:
            small_output_f = False
        if args["verbosity"] == '0':
            df["contigA"] = [''] * len(df)
            df["contigB"] = [''] * len(df)
        elif args["verbosity"] == '1':
            has_alt = [True if isinstance(i, str) and i[0] != '<' else False for i in df['variant_seq']]
            df["contigA"] = ['' if a else c for c, a in zip(df['contigA'], has_alt)]
            df["contigB"] = ['' if a else c for c, a in zip(df['contigB'], has_alt)]

        n_fields = len(col_names(small_output_f)[-1])
    else:
        small_output_f = False

    for idx, r in df.iterrows():
        if idx in seen_idx:
            continue

        format_f, df_rows = gen_format_fields(r, df, names, n_fields, small_output_f)
        if "partners" in r and r["partners"] is not None and r["partners"] != ".":
            seen_idx |= set(r["partners"])

        r_main = make_main_record(r, dysgu_version, count, format_f, df_rows, add_kind, small_output_f)
        recs.append(r_main)
        count += 1

    if sort_output:
        for rec in sorted(recs, key=lambda x: (x[0], x[1])):
            outfile.write("\t".join(list(map(str, rec))) + "\n")
    else:
        for rec in recs:
            outfile.write("\t".join(list(map(str, rec))) + "\n")

    return count
