#!python
#cython: language_level=2, boundscheck=True, wraparound=True
#distutils: language=c++

import numpy as np
cimport numpy as np
import logging
from map_set_utils import echo
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
        logging.info("Reading {}:{}-{}".format(c, s, e))
        for r in bam.fetch(c, int(s), int(e)):
            yield r


cpdef list col_names(extended, small_output):  # todo fix for no-contigs view command

    if small_output:
        return ["chrA", "posA", "chrB", "posB", "sample", "id", "grp_id", "n_in_grp", "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B", "svlen", "svlen_precise", "rep", "gc",
          ["GT", "GQ", "MAPQpri", "su", "spanning", "pe", "supp", "sc", "bnd",
         "raw_reads_10kb", "neigh10kb", "plus",
                "minus", "remap_score", "remap_ed", "bad_clip_count", "fcc", "inner_cn", "outer_cn", "prob"]
            ]
    if extended:  # to do fix this
        return ["chrA", "posA", "chrB", "posB", "sample", "id", "grp_id", "n_in_grp",  "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B", 'contigA', 'contigB', "svlen", "svlen_precise",  "rep", "gc",
         ["GT", "GQ", "DP", "DN", "DApri", "DAsupp",  "NMpri", "NMsupp", "NMbase", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "su", "spanning", "pe", "supp", "sc", "bnd", "sqc", "block_edge",
         "raw_reads_10kb", "mcov",
          "linked", "neigh", "neigh10kb",  "ref_bases", "plus",
                "minus", "n_gaps", "n_sa", "n_xa", "n_unmapped_mates", "double_clips", "remap_score", "remap_ed", "bad_clip_count", "fcc", "n_small_tlen", "ras", "fas",
                "inner_cn", "outer_cn", "compress", "ref_rep", "prob"]
            ]
    else:
        return ["chrA", "posA", "chrB", "posB", "sample", "id", "grp_id", "n_in_grp", "kind", "type", "svtype", "join_type", "cipos95A", "cipos95B", 'contigA', 'contigB', "svlen", "svlen_precise",  "rep", "gc",
          ["GT", "GQ", "NMpri", "NMsupp", "NMbase", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "su", "spanning", "pe", "supp", "sc", "bnd", "sqc", "block_edge",
         "raw_reads_10kb", "mcov",
          "linked", "neigh", "neigh10kb",  "ref_bases", "plus",
                "minus", "n_gaps", "n_sa", "n_xa", "n_unmapped_mates", "double_clips", "remap_score", "remap_ed", "bad_clip_count", "fcc", "n_small_tlen", "ras", "fas",
                "inner_cn", "outer_cn", "compress", "ref_rep", "prob"]
            ]


def make_main_record(r, version, index, format_f, df_rows, add_kind, extended, small_output):

    # Pick best row (best support, or highest prob if available
    rep, repsc, lenprec = 0, 0, 1
    mean_prob, max_prob = None, None
    if len(format_f) > 1:

        best = sorted([(int(v["su"]), k) for k, v in df_rows.items()], reverse=True)[0][1]
        probs = [v["prob"] for k, v in df_rows.items()]
        mean_prob = np.mean(probs)
        max_prob = np.max(probs)
        r = df_rows[best]
        gc = round(r["gc"], 2)
        if not small_output:
            rep = r["rep"]
            repsc = r["rep_sc"]
            lenprec = 1 if "svlen_precise" not in r else r["svlen_precise"]
        n_expansion = r["n_expansion"]
        stride = r["stride"]
        exp_seq = r["exp_seq"]
        ref_poly = r["ref_poly_bases"]
        # q_gaps = r["query_gap"]
        overlaps = r["query_overlap"]

        su, pe, sr, sc, bnd, wr = 0, 0, 0, 0, 0, 0
        # probs = []
        for row in df_rows.values():
            pe += row["pe"]
            sr += row["supp"]
            sc += row["sc"]
            su += row["su"]
            bnd += row["bnd"]
            wr += row["spanning"]
            # probs.append(row["Prob"])
        # probs = round(np.median(probs), 3)

    else:
        pe = r["pe"]
        sr = r["supp"]
        sc = r["sc"]
        su = r["su"]
        bnd = r["bnd"]
        wr = r["spanning"]
        # probs = r["Prob"]
        gc = round(r["gc"], 2)
        if not small_output:
            rep = r["rep"]
            repsc = r["rep_sc"]
            lenprec = 1 if "svlen_precise" not in r else r["svlen_precise"]
        n_expansion = r["n_expansion"]
        stride = r["stride"]
        exp_seq = r["exp_seq"]
        ref_poly = r["ref_poly_bases"]
        # q_gaps = r["query_gap"]
        overlaps = r["query_overlap"]

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

    if not small_output:
        info_extras += [
                    f"REP={'%.3f' % float(rep)}",
                    f"REPSC={'%.3f' % float(repsc)}",

                    ]

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

    if mean_prob is not None:
        info_extras += [f"MeanPROB={round(mean_prob, 3)}", f"MaxPROB={round(max_prob, 3)}"]

    if small_output:
        fmt_keys = "GT:GQ:MAPQP:SU:WR:PE:SR:SC:BND:COV:NEIGH10:PS:MS:RMS:RED:BCC:FCC:ICN:OCN:PROB"
        # if "prob" in r:
        #     fmt_keys += ":PROB"

    elif extended:
        fmt_keys = "GT:GQ:DP:DN:DAP:DAS:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BND:SQC:SCW:BE:COV:MCOV:LNK:NEIGH:NEIGH10:RB:PS:MS:NG:NSA:NXA:NMU:NDC:RMS:RED:BCC:FCC:STL:RAS:FAS:ICN:OCN:CMP:RR:JIT:PROB"
        # if "prob" in r:
        #     fmt_keys += ":PROB"
    else:
        fmt_keys = "GT:GQ:NMP:NMS:NMB:MAPQP:MAPQS:NP:MAS:SU:WR:PE:SR:SC:BND:SQC:SCW:BE:COV:MCOV:LNK:NEIGH:NEIGH10:RB:PS:MS:NG:NSA:NXA:NMU:NDC:RMS:RED:BCC:FCC:STL:RAS:FAS:ICN:OCN:CMP:RR:JIT:PROB"
        # if "prob" in r:
        #     fmt_keys += ":PROB"

    rec = [r["chrA"], r["posA"], index, ".", f"<{r['svtype']}>", ".", "." if "filter" not in r else r['filter'],
           # INFO line
           ";".join([f"SVMETHOD=DYSGUv{version}",
                   f"SVTYPE={r['svtype']}",
                   f"END={r['posB']}",
                   f"CHR2={r['chrB']}",
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

    return rec


def get_fmt(r, extended, small_output):
    if small_output:
        v = [r["GT"], r["GQ"], r['MAPQpri'],
                                  r['su'], r['spanning'], r['pe'], r['supp'],
                                  r['sc'], r['bnd'], r['raw_reads_10kb'], r['neigh10kb'],
                                  r["plus"], r["minus"], r["remap_score"], r["remap_ed"], r["bad_clip_count"], round(r["fcc"], 3),
                                round(r["inner_cn"], 3), round(r["outer_cn"], 3), r['prob']
             ]
        return v

    elif extended:
        v = [r["GT"], r['GQ'], r['DP'], r['DN'], r['DApri'], r['DAsupp'], r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
                                      r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
                                      r['sc'], r['bnd'], round(r['sqc'], 2), round(r['scw'], 1), r['block_edge'], r['raw_reads_10kb'], round(r['mcov'], 2), int(r['linked']), r['neigh'], r['neigh10kb'],
                                      r['ref_bases'], r["plus"], r["minus"], r['n_gaps'], round(r["n_sa"], 2), round(r["n_xa"], 2),
                                      round(r["n_unmapped_mates"], 2), r["double_clips"], r["remap_score"], r["remap_ed"], r["bad_clip_count"], round(r["fcc"], 3), r["n_small_tlen"], r["ras"], r['fas'],
                                    round(r["inner_cn"], 3), round(r["outer_cn"], 3), round(r["compress"], 2), round(r["ref_rep"], 3), round(r["jitter"], 3), r['prob']

             ]
        return v

    else:
        v = [r["GT"], r["GQ"], r['NMpri'], r['NMsupp'], r['NMbase'], r['MAPQpri'],
                                  r['MAPQsupp'], r['NP'], r['maxASsupp'], r['su'], r['spanning'], r['pe'], r['supp'],
                                  r['sc'], r['bnd'], round(r['sqc'], 2), round(r['scw'], 1), r['block_edge'], r['raw_reads_10kb'], round(r['mcov'], 2), int(r['linked']), r['neigh'], r['neigh10kb'],
                                  r['ref_bases'], r["plus"], r["minus"], r['n_gaps'], round(r["n_sa"], 2),
                                  round(r["n_xa"], 2), round(r["n_unmapped_mates"], 2), r["double_clips"], r["remap_score"], r["remap_ed"], r["bad_clip_count"], round(r["fcc"], 3), r["n_small_tlen"], r["ras"], r['fas'],
                                round(r["inner_cn"], 3), round(r["outer_cn"], 3), round(r["compress"], 2), round(r["ref_rep"], 3), round(r["jitter"], 3), r['prob']
             ]
        return v


def gen_format_fields(r, df, names, extended, n_fields, small_output):

    if len(names) == 1:
        return {0: get_fmt(r, extended, small_output)}, {}

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
            format_fields[name] = get_fmt(row, extended, small_output)
        else:
            format_fields[name] = [0] * n_fields

    return format_fields, cols



def to_vcf(df, args, names, outfile, show_names=True,  contig_names="", extended_tags=False, header=None,
           small_output_f=True):
    if header is None:
        if extended_tags:
            HEADER = """##fileformat=VCFv4.2
##source=DYSGU
##FILTER=<ID=lowProb,Description="Probability below threshold set with --thresholds">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
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
##INFO=<ID=MeanPROB,Number=1,Type=Float,Description="Mean probability of event being true across samples">
##INFO=<ID=MaxPROB,Number=1,Type=Float,Description="Max probability of event being true across samples">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality phred scaled">
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
##FORMAT=<ID=BND,Number=1,Type=Integer,Description="Number of break-end alignments supporting the variant">
##FORMAT=<ID=SQC,Number=1,Type=Float,Description="Soft-clip quality value correlation between reads">
##FORMAT=<ID=SCQ,Number=1,Type=Float,Description="Soft-clip quality weight value">
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Mean read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=MCOV,Number=1,Type=Float,Description="Maximum read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other break points within 1 bp of break site">
##FORMAT=<ID=NEIGH10,Number=1,Type=Integer,Description="Number of other break points within 10 kp of break site">
##FORMAT=<ID=RB,Number=1,Type=Integer,Description="Number of reference bases in contigs">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Number of reads on plus strand">
##FORMAT=<ID=MS,Number=1,Type=Integer,Description="Number of reads on minus strand">
##FORMAT=<ID=NG,Number=1,Type=Float,Description="Mean number of small gaps < 30 bp">
##FORMAT=<ID=NSA,Number=1,Type=Float,Description="Mean number of SA tags per read">
##FORMAT=<ID=NXA,Number=1,Type=Float,Description="Mean number of XA tags per read">
##FORMAT=<ID=NMU,Number=1,Type=Float,Description="Mean number of mates unmapped per read">
##FORMAT=<ID=NDC,Number=1,Type=Integer,Description="Number of double-clips, alignments with left and right clips">
##FORMAT=<ID=RMS,Number=1,Type=Integer,Description="Remapping score">
##FORMAT=<ID=RED,Number=1,Type=Integer,Description="Remapping edit distance">
##FORMAT=<ID=BCC,Number=1,Type=Integer,Description="Bad soft-clip count within +/- 500 bp">
##FORMAT=<ID=FCC,Number=1,Type=Float,Description="Fold-coverage change for SVs > 200 bp">
##FORMAT=<ID=STL,Number=1,Type=Integer,Description="N reads with small TLEN below 0.05% of distribution">
##FORMAT=<ID=RAS,Number=1,Type=Integer,Description="Reverse soft-clip to alignment score">
##FORMAT=<ID=FAS,Number=1,Type=Integer,Description="Forward soft-clip to alignment score">
##FORMAT=<ID=ICN,Number=1,Type=Float,Description="Inner copy number">
##FORMAT=<ID=OCN,Number=1,Type=Float,Description="Outer copy number">
##FORMAT=<ID=CMP,Number=1,Type=Float,Description="Compression ratio of contigs">
##FORMAT=<ID=RR,Number=1,Type=Float,Description="Repeat score for reference">
##FORMAT=<ID=JIT,Number=1,Type=Float,Description="SV length jitter">
##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event being true">{}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

        else:
            HEADER = """##fileformat=VCFv4.2
##source=DYSGU
##FILTER=<ID=lowProb,Description="Probability below threshold set with --thresholds">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
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
##INFO=<ID=MeanPROB,Number=1,Type=Float,Description="Mean probability of event being true across samples">
##INFO=<ID=MaxPROB,Number=1,Type=Float,Description="Max probability of event being true across samples">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality phred scaled">
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
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Mean read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=MCOV,Number=1,Type=Float,Description="Maximum read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other break points within 1 bp of break site">
##FORMAT=<ID=NEIGH10,Number=1,Type=Integer,Description="Number of other break points within 10 kp of break site">
##FORMAT=<ID=RB,Number=1,Type=Integer,Description="Number of reference bases in contigs">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Number of reads on plus strand">
##FORMAT=<ID=MS,Number=1,Type=Integer,Description="Number of reads on minus strand">
##FORMAT=<ID=NG,Number=1,Type=Float,Description="Mean number of small gaps < 30 bp">
##FORMAT=<ID=NSA,Number=1,Type=Float,Description="Mean number of SA tags per read">
##FORMAT=<ID=NXA,Number=1,Type=Float,Description="Mean number of XA tags per read">
##FORMAT=<ID=NMU,Number=1,Type=Float,Description="Mean number of mates unmapped per read">
##FORMAT=<ID=NDC,Number=1,Type=Integer,Description="Number of double-clips, alignments with left and right clips">
##FORMAT=<ID=RMS,Number=1,Type=Integer,Description="Remapping score">
##FORMAT=<ID=RED,Number=1,Type=Integer,Description="Remapping edit distance">
##FORMAT=<ID=BCC,Number=1,Type=Integer,Description="Bad soft-clip count within +/- 500 bp">
##FORMAT=<ID=FCC,Number=1,Type=Float,Description="Fold-coverage change for SVs > 200 bp">
##FORMAT=<ID=STL,Number=1,Type=Integer,Description="N reads with small TLEN below 0.05% of distribution">
##FORMAT=<ID=RAS,Number=1,Type=Integer,Description="Reverse soft-clip to alignment score">
##FORMAT=<ID=FAS,Number=1,Type=Integer,Description="Forward soft-clip to alignment score">
##FORMAT=<ID=ICN,Number=1,Type=Float,Description="Inner copy number">
##FORMAT=<ID=OCN,Number=1,Type=Float,Description="Outer copy number">
##FORMAT=<ID=CMP,Number=1,Type=Float,Description="Compression ratio of contigs">
##FORMAT=<ID=RR,Number=1,Type=Float,Description="Repeat score for reference">
##FORMAT=<ID=JIT,Number=1,Type=Float,Description="SV length jitter">
##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event being true">{}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

# ##INFO=<ID=MPROB,Number=1,Type=Float,Description="Median probability of event across samples">
# ##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event">

    else:
        HEADER = header

    outfile.write(HEADER.format(contig_names) + "\t" + "\t".join(names) + "\n")

    if show_names:
        logging.info("Input samples: {}".format(str(list(names))))

    version = pkg_resources.require("dysgu")[0].version

    seen_idx = set([])
    cnames = ['raw_reads_10kb', 'NMpri', 'NMsupp', 'MAPQpri', 'MAPQsupp', "NMbase", "n_gaps"]
    if extended_tags:
        cnames += ['DP', 'DN', 'DApri', 'DAsupp']

    for col in cnames:
        if col in df.columns:
            df[col] = df[col].astype(float).round(2)

    for col in ['maxASsupp', 'neigh', 'neigh10kb']:
        if col in df.columns:
            df[col] = [int(i) if i == i else 0 for i in df[col]]

    count = 0
    recs = []
    jobs = []

    add_kind = args["add_kind"] == "True"
    if args["metrics"]:
        small_output_f = False

    n_fields = len(col_names(extended_tags, small_output_f)[-1])
    for idx, r in df.iterrows():

        if idx in seen_idx:
            continue

        format_f, df_rows = gen_format_fields(r, df, names, extended_tags, n_fields, small_output_f)
        if "partners" in r and r["partners"] is not None:
            seen_idx |= set(r["partners"])

        r_main = make_main_record(r, version, count, format_f, df_rows, add_kind, extended_tags, small_output_f)

        recs.append(r_main)
        count += 1

    for rec in sorted(recs, key=lambda x: (x[0], x[1])):
        outfile.write("\t".join(list(map(str, rec))) + "\n")
    return count
