#!python
#cython: language_level=2, boundscheck=False, wraparound=False
#distutils: language=c++

import numpy as np
cimport numpy as np
cimport cython
import click
from collections import defaultdict
import ncls
import pkg_resources
import sortedcontainers
import pandas as pd
import os


DTYPE = np.float
ctypedef np.float_t DTYPE_t


from libc.stdlib cimport malloc
import re

cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0',
                       '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0',  'a',  'a' ]

np.random.seed(0)


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


cdef list get_bed_regions(str bed):
    b = [tuple([int(j) if j.isdigit() else j for j in i.strip().split("\t")[:3]]) for i in open(bed, "r")
         if i[0] != "#" and len(i) > 0 and "\t" in i]
    if len(b) == 0:
        raise ValueError("Bed regions not formatted correctly")
    return b


cpdef dict overlap_regions(str bed, int_chroms=False, infile=None):
    if not bed:
        return {}
    regions = get_bed_regions(bed)
    chrom_interval_start = defaultdict(list)
    chrom_interval_end = defaultdict(list)
    for c, s, e in regions:
        if int_chroms:
            c = infile.gettid(c)
        chrom_interval_start[c].append(int(s))
        chrom_interval_end[c].append(int(e))

    regions = {k: ncls.NCLS(np.array(chrom_interval_start[k]),
                            np.array(chrom_interval_end[k]),
                            np.array(chrom_interval_start[k])) for k in chrom_interval_start}

    return regions


cpdef int intersecter_int_chrom(dict tree, int chrom, int start, int end):
    if not tree:
        return 0
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return 1
        else:
            return 0
    else:
        return 0


cpdef int intersecter_str_chrom(dict tree, str chrom, int start, int end):
    if not tree:
        return 0
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return 1
        else:
            return 0
    else:
        return 0


def mk_dest(d):
    if d is not None and not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            raise OSError("Couldn't create directory {}".format(d))



# def get_bed_regions(bed):
#     b = [tuple([int(j) if j.isdigit() else j for j in i.strip().split("\t")[:3]]) for i in open(bed, "r")
#          if i[0] != "#" and len(i) > 0 and "\t" in i]
#     if len(b) == 0:
#         raise ValueError("Bed regions not formatted correctly")
#     return b


def overlap_regionspy(bed):
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


def intersecterpy(tree, chrom, start, end):
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


cpdef list col_names(extended):
    if extended:
        return ["chrA", "posA", "chrB", "posB", "sample", "id", "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B",
         "DP", "DN", "DApri", "DAsupp",  "NMpri", "NMsupp", "NMbase", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "su", "pe", "supp", "sc", "block_edge",
         "raw_reads_10kb",
          "linked", "contigA", "contigB",  "gc", "neigh", "neigh10kb", "rep", "rep_sc", "svlen_precise", "ref_bases", "svlen", "plus",
                "minus", "n_gaps", "n_sa", "n_xa", "n_unmapped_mates", "double_clips", "switched"
            ]
    else:
        return ["chrA", "posA", "chrB", "posB", "sample", "id", "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B",
          "NMpri", "NMsupp", "NMbase", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "su", "pe", "supp", "sc", "block_edge",
         "raw_reads_10kb",
          "linked", "contigA", "contigB",  "gc", "neigh", "neigh10kb", "rep", "rep_sc", "svlen_precise", "ref_bases", "svlen", "plus",
                "minus", "n_gaps", "n_sa", "n_xa", "n_unmapped_mates", "double_clips", "switched"
            ]


def make_main_record(r, version, index, format_f, df_rows, add_kind, extended):

    # Pick best row (best support, or highest prob if available
    if len(format_f) > 1:

        best = sorted([(int(v["su"]), k) for k, v in df_rows.items()], reverse=True)[0][1]
        r = df_rows[best]
        gc = r["gc"]
        rep = r["rep"]
        repsc = r["rep_sc"]
        lenprec = 1 if "svlen_precise" not in r else r["svlen_precise"]
        n_expansion = r["n_expansion"]
        stride = r["stride"]
        exp_seq = r["exp_seq"]
        ref_poly = r["ref_poly_bases"]

        su, pe, sr, sc, wr = 0, 0, 0, 0, 0
        # probs = []
        for row in df_rows.values():
            pe += row["pe"]
            sr += row["supp"]
            sc += row["sc"]
            su += row["su"]
            wr += row["spanning"]
            # probs.append(row["Prob"])
        # probs = round(np.median(probs), 3)

    else:
        pe = r["pe"]
        sr = r["supp"]
        sc = r["sc"]
        su = r["su"]
        wr = r["spanning"]
        # probs = r["Prob"]
        gc = r["gc"]
        rep = r["rep"]
        repsc = r["rep_sc"]
        lenprec = 1 if "svlen_precise" not in r else r["svlen_precise"]
        n_expansion = r["n_expansion"]
        stride = r["stride"]
        exp_seq = r["exp_seq"]
        ref_poly = r["ref_poly_bases"]

    samp = r["sample"]
    read_kind = r["type"]

    if r["chrA"] == r["chrB"] and int(r["posA"]) > int(r["posB"]):
        chrA, posA, cipos95A, contig2 = r["chrA"], r["posA"], r["cipos95A"], r["contigB"]
        r["chrA"] = r["chrB"]
        r["posA"] = r["posB"]
        r["cipos95A"] = r["cipos95B"]
        r["chrB"] = chrA
        r["posB"] = posA
        r["cipos95B"] = cipos95A
        r["contigB"] = r["contigA"]
        r["contigA"] = contig2

    info_extras = []
    if r["chrA"] == r["chrB"]:
        # svlen = abs(r["posA"] - r["posB"])
        if 'svlen' in r:
            info_extras.append(f"SVLEN={r['svlen']}")
        else:
            info_extras.append(f"SVLEN=0")

    if r["contigA"]:
        info_extras.append(f"CONTIGA={r['contigA']}")
    if r["contigB"]:
        info_extras.append(f"CONTIGB={r['contigB']}")

    if add_kind:
        info_extras += [f"KIND={r['kind']}"]

    info_extras += [f"GC={gc}",
                    f"REP={'%.3f' % float(rep)}",
                    f"REPSC={'%.3f' % float(repsc)}",
                    f"LPREC={lenprec}",
                    f"NEXP={n_expansion}",
                    f"STRIDE={stride}",
                    f"EXPSEQ={exp_seq}",
                    f"RPOLY={ref_poly}",
                    f"SU={su}",
                    f"WR={wr}",
                    f"PE={pe}",
                    f"SR={sr}",
                    f"SC={sc}",
                    f"RT={read_kind}"]

    if extended:
        fmt_keys = "GT:DP:DN:DAP:DAS:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BE:COV:LNK:NEIGH:NEIGH10:RB:PS:MS:NG:NSA:NXA:NMU:NDC:SW"
        if "prob" in r:
            fmt_keys += ":PROB"
    else:
        fmt_keys = "GT:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BE:COV:LNK:NEIGH:NEIGH10:RB:PS:MS:NG:NSA:NXA:NMU:NDC:SW"
        if "prob" in r:
            fmt_keys += ":PROB"

    rec = [r["chrA"], r["posA"], index, ".", f"<{r['svtype']}>", ".", ".",
           # INFO line
           ";".join([f"SVMETHOD=DYSGUv{version}",
                   f"SVTYPE={r['svtype']}",
                   f"END={r['posB']}",
                   f"CHR2={r['chrB']}",
                   f"CT={r['join_type']}",
                   f"CIPOS95={r['cipos95A']}",
                   f"CIEND95={r['cipos95B']}",

                   ] + info_extras),
           fmt_keys  # :PROB
           ]
    # FORMAT line(s)
    for item in format_f.values():
        rec.append(":".join(map(str, item)))

    return rec


def get_fmt(r, extended):
    if extended:

        v = ["./.", r['DP'], r['DN'], r['DApri'], r['DAsupp'], r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
                                      r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
                                      r['sc'], r['block_edge'], r['raw_reads_10kb'], r['linked'], r['neigh'], r['neigh10kb'],
                                      r['ref_bases'], r["plus"], r["minus"], r['n_gaps'], round(r["n_sa"], 2), round(r["n_xa"], 2),
                                      round(r["n_unmapped_mates"], 2), r["double_clips"], r["switched"]]
        if "prob" in r:
            v.append(r["prob"])
        return v

    else:
        v = ["./.", r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
                                  r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
                                  r['sc'], r['block_edge'], r['raw_reads_10kb'], r['linked'], r['neigh'], r['neigh10kb'],
                                  r['ref_bases'], r["plus"], r["minus"], r['n_gaps'], round(r["n_sa"], 2),
                                  round(r["n_xa"], 2), round(r["n_unmapped_mates"], 2), r["double_clips"], r["switched"]]
        if "prob" in r:
            v.append(r["prob"])
        return v


def gen_format_fields(r, df, names, extended, n_fields):

    if len(names) == 1:
        return {0: get_fmt(r, extended)}, {}

    cols = {}
    if "partners" in r:
        if not isinstance(r["partners"], set):
            if len(r["partners"]) == 0 or pd.isna(r["partners"]):
                r["partners"] = []
            else:
                r["partners"] = [int(i.split(",")[1]) for i in r["partners"].split("|")]

        for idx in r["partners"]:
            r2 = df.loc[idx]  # iloc
            cols[r2["table_name"]] = r2

    if "table_name" in r:
        cols[r["table_name"]] = r

    format_fields = sortedcontainers.SortedDict()
    for name in names:
        if name in cols:
            row = cols[name]
            format_fields[name] = get_fmt(row, extended)
        else:
            format_fields[name] = [0] * n_fields

    return format_fields, cols



def to_vcf(df, args, names, outfile, n_fields=19, show_names=True,  contig_names="", extended_tags=False, header=None):
    if header is None:
        if extended_tags:
            HEADER = """##fileformat=VCFv4.2
##source=DYSGU
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
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
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean distance-to-pair metric supporting the variant">
##FORMAT=<ID=DN,Number=1,Type=Float,Description="Mean distance-to-normal metric supporting the variant">
##FORMAT=<ID=DAP,Number=1,Type=Float,Description="Mean distance-to-alignment metric for primary alignments">
##FORMAT=<ID=DAS,Number=1,Type=Float,Description="Mean distance-to-alignment metric for supplementary alignments">
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
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Maximum read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other break points within 1 bp of break site">
##FORMAT=<ID=NEIGH10,Number=1,Type=Integer,Description="Number of other break points within 10 kp of break site">
##FORMAT=<ID=RB,Number=1,Type=Integer,Description="Number of reference bases in contigs">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Number of reads on plus strand">
##FORMAT=<ID=MS,Number=1,Type=Integer,Description="Number of reads on minus strand">
##FORMAT=<ID=NG,Number=1,Type=Integer,Description="Mean number of small gaps < 30 bp">
##FORMAT=<ID=NSA,Number=1,Type=Integer,Description="Mean number of SA tags per read">
##FORMAT=<ID=NXA,Number=1,Type=Integer,Description="Mean number of XA tags per read">
##FORMAT=<ID=NMU,Number=1,Type=Integer,Description="Mean number of mates unmapped per read">
##FORMAT=<ID=NDC,Number=1,Type=Integer,Description="Number of double-clips, alignments with left and right clips">
##FORMAT=<ID=SW,Number=1,Type=Integer,Description="Sv-type was switched by remapping">{}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

        else:
            HEADER = """##fileformat=VCFv4.2
##source=DYSGU
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
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
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
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
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Maximum read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other break points within 1 bp of break site">
##FORMAT=<ID=NEIGH10,Number=1,Type=Integer,Description="Number of other break points within 10 kp of break site">
##FORMAT=<ID=RB,Number=1,Type=Integer,Description="Number of reference bases in contigs">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Number of reads on plus strand">
##FORMAT=<ID=MS,Number=1,Type=Integer,Description="Number of reads on minus strand">
##FORMAT=<ID=NG,Number=1,Type=Integer,Description="Mean number of small gaps < 30 bp">
##FORMAT=<ID=NSA,Number=1,Type=Integer,Description="Mean number of SA tags per read">
##FORMAT=<ID=NXA,Number=1,Type=Integer,Description="Mean number of XA tags per read">
##FORMAT=<ID=NMU,Number=1,Type=Integer,Description="Mean number of mates unmapped per read">
##FORMAT=<ID=NDC,Number=1,Type=Integer,Description="Number of double-clips, alignments with left and right clips">
##FORMAT=<ID=SW,Number=1,Type=Integer,Description="Sv-type was switched by remapping">{}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

# ##INFO=<ID=MPROB,Number=1,Type=Float,Description="Median probability of event across samples">
# ##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event">

    else:
        HEADER = header

    outfile.write(HEADER.format(contig_names) + "\t" + "\t".join(names) + "\n")

    if show_names:
        click.echo("Input samples: {}".format(str(list(names))), err=True)

    version = pkg_resources.require("dysgu")[0].version

    # if len(names) > 1:
    #     dm = df.sort_values(["partners"], ascending=False)
    # else:
    #     dm = df

    seen_idx = set([])

    cnames = ['raw_reads_10kb', 'NMpri', 'NMsupp', 'MAPQpri', 'MAPQsupp', "NMbase", "n_gaps"]
    if extended_tags:
        # cnames = ['raw_reads_10kb',  'NMpri', 'NMsupp', 'MAPQpri', 'MAPQsupp']
        cnames += ['DP', 'DN', 'DApri', 'DAsupp',]

    for col in cnames:
        df[col] = df[col].astype(float).round(2)

    for col in ['maxASsupp', 'neigh', 'neigh10kb']:
        df[col] = [int(i) for i in df[col]]

    count = 0
    recs = []
    jobs = []

    add_kind = args["add_kind"] == "True"
    for idx, r in df.iterrows():

        if idx in seen_idx:
            continue

        format_f, df_rows = gen_format_fields(r, df, names, extended_tags, n_fields)
        # click.echo(format_f, err=True)
        if "partners" in r:
            seen_idx |= set(r["partners"])

        r_main = make_main_record(r, version, count, format_f, df_rows, add_kind, extended_tags)
        # click.echo(r_main, err=True)
        # quit()
        recs.append(r_main)
        count += 1

    for rec in sorted(recs, key=lambda x: (x[0], x[1])):
        outfile.write("\t".join(list(map(str, rec))) + "\n")
    return count
