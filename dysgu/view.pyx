#cython: language_level = 3

import os
import pandas as pd
import click
import sys
import sortedcontainers
from dysgu import cluster
import pkg_resources
import numpy as np
import time
import datetime
from collections import defaultdict


HEADER = """##fileformat=VCFv4.2
##source=DYSGU
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=CIPOS95,Number=1,Type=Integer,Description="Confidence interval size (95%) around POS for imprecise variants">
##INFO=<ID=CIEND95,Number=1,Type=Integer,Description="Confidence interval size (95%) around END for imprecise variants">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=KIND,Number=1,Type=String,Description="Kind of join with respect to input regions">
##INFO=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of supplementary reads supporting the variant across all samples">
##INFO=<ID=SC,Number=1,Type=Integer,Description="Number of soft-clip reads supporting the variant across all samples">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Contig from CHROM POS">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Contig from CHR2 END">
##INFO=<ID=GC,Number=1,Type=Float,Description="GC% of assembled contigs">
##INFO=<ID=REP,Number=1,Type=Float,Description="Homopolymer repeat score for contigs">
##INFO=<ID=RB,Number=1,Type=Float,Description="Reference bases assembled">
##INFO=<ID=MPROB,Number=1,Type=Float,Description="Median probability of event across samples">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean distance-to-pair metric supporting the variant">
##FORMAT=<ID=DN,Number=1,Type=Float,Description="Mean distance-to-normal metric supporting the variant">
##FORMAT=<ID=DAP,Number=1,Type=Float,Description="Mean distance-to-alignment metric for primary alignments">
##FORMAT=<ID=DAS,Number=1,Type=Float,Description="Mean distance-to-alignment metric for supplementary alignments">
##FORMAT=<ID=NMP,Number=1,Type=Float,Description="Mean edit distance for primary alignments supporting the variant">
##FORMAT=<ID=NMS,Number=1,Type=Float,Description="Mean edit distance for supplementary alignments supporting the variant">
##FORMAT=<ID=MAPQP,Number=1,Type=Float,Description="Mean MAPQ for primary reads supporting the variant">
##FORMAT=<ID=MAPQS,Number=1,Type=Float,Description="Mean MAPQ for supplementary reads supporting the variant">
##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Number of alignments in normal-pair orientation supporting the variant">
##FORMAT=<ID=MAS,Number=1,Type=Integer,Description="Maximum alignment score of supplementary reads supporting the variant">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supplementary alignments supporting the variant">
##FORMAT=<ID=SC,Number=1,Type=Integer,Description="Number of soft-clipped alignments supporting the variant">
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Mean read coverage +/- 10kb around break sites">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other beak points within 100 bp or break sites">
##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""


def echo(*args):
    click.echo(args, err=True)


def open_outfile(args):

    if args["separate"] == "True":
        outfiles = {}
        for item in args["input_files"]:
            name, ext = os.path.splitext(item)
            outname = f"{name}.{args['post_fix']}.csv"
            outfiles[name] = open(outname, "w")
        return outfiles

    if args["svs_out"] == "-" or args["svs_out"] is None:
        click.echo("SVs output to stdout", err=True)
        outfile = sys.stdout
    else:
        click.echo("SVs output to {}".format(args["svs_out"]), err=True)
        outfile = open(args["svs_out"], "w")

    return outfile


def merge_df(df, tree=None, merge_within_sample=False):

    df.reset_index(inplace=True)
    potential = []
    for idx, r in df.iterrows():
        r = dict(r)
        r.update({"event_id": idx, "contig": r["contigA"], "contig2": r["contigB"],})
        potential.append(r)

    if not merge_within_sample:
        found = cluster.merge_events(potential, 15, tree, try_rev=False, pick_best=True, add_partners=True)
        ff = defaultdict(list)

        for f in found:

            if "partners" not in f:
                ff[f["event_id"]] = []
            else:
                # Remove partners from same sample
                current = f["table_name"]

                targets = set([])
                # others = []
                for item in f["partners"]:  # Only merge with one row per sample
                    t_name = df.iloc[item]["table_name"]
                    if t_name != current and t_name not in targets:
                        ff[f["event_id"]].append(item)
                        ff[item].append(f["event_id"])

                        targets.add(t_name)
                # echo(f["event_id"], f["partners"], others, targets)
                # echo(f["event_id"], f["partners"], targets)
                #ff[f["event_id"]] = others
        df["partners"] = [ff[i] if i in ff else [] for i in df.index]
        return df
    else:
        found = cluster.merge_events(potential, 15, tree, try_rev=False, pick_best=True, add_partners=False)
        return pd.DataFrame.from_records(found)




def make_main_record(r, version, index, format_f, df_rows, df):

    # Pick best row (best support, or highest prob if available
    if len(format_f) > 1:
        best = sorted([(v["Prob"], v["pe"] + v["supp"], k) for k, v in df_rows.items()],
                      key=lambda x: (-x[0], -x[1]))[0][2]
        r = df_rows[best]

        su, pe, sr, sc = 0, 0, 0, 0
        probs = []
        for row in df_rows.values():
            pe += row["pe"]
            sr += row["supp"]
            sc += row["sc"]
            su += (row["pe"] + row["supp"])
            probs.append(row["Prob"])
        probs = round(np.median(probs), 3)

    else:
        pe = r["pe"]
        sr = r["supp"]
        sc = r["sc"]
        su = (r["pe"] + r["supp"])
        probs = r["Prob"]

    samp = r["sample"]

    info_extras = []
    if r["chrA"] == r["chrB"]:
        svlen = abs(r["posA"] - r["posB"])
        info_extras.append(f"SVLEN={svlen}")

    if r["contigA"]:
        info_extras.append(f"CONTIGA={r['contigA']}")
    if r["contigB"]:
        info_extras.append(f"CONTIGB={r['contigB']}")

    info_extras += [f"SU={su}",
                    f"PE={pe}",
                    f"SR={sr}",
                    f"SC={sc}",
                    f"MPROB={probs}"]

    # if "partners" in r:
    #     p = "|".join([f"{df.iloc[i]['table_name']},{df.iloc[i]['id']}" for i in r["partners"]])
    #     if p:
    #         info_extras += [f"PARTNERS={p}"]

    rec = [r["chrA"], r["posA"], index, ".", f"<{r['svtype']}>", ".", ".",
           # INFO line
           ";".join([f"SVMETHOD=DYSGUv{version}",
                   f"SVTYPE={r['svtype']}",
                   f"END={r['posB']}",
                   f"CHR2={r['chrB']}",
                   f"CT={r['join_type']}",
                   f"CIPOS95={r['cipos95A']}",
                   f"CIEND95={r['cipos95B']}",
                   f"KIND={r['kind']}",

                   ] + info_extras),
           "DP:DN:DAP:DAS:NMP:NMS:MAPQP:MAPQS:NP:MAS:SU:PE:SR:SC:BE:COV:LNK:NEIGH:RB:PROB"
           ]
    # FORMAT line(s)
    for item in format_f.values():
        rec.append(":".join(map(str, item)))


    return rec


def gen_format_fields(r, df, names):

    cols = {}
    if "partners" in r:
        if not isinstance(r["partners"], list):
            if len(r["partners"]) == 0 or pd.isna(r["partners"]):
                r["partners"] = []
            else:
                r["partners"] = [int(i.split(",")[1]) for i in r["partners"].split("|")]
        # echo(r, r["partners"])
        for idx in r["partners"]:

            r2 = df.iloc[idx]
            cols[r2["table_name"]] = r2

    cols[r["table_name"]] = r

    format_fields = sortedcontainers.SortedDict()

    for name in names:

        if name in cols:
            format_fields[name] = ([r['DP'], r['DN'], r['DApri'], r['DAsupp'], r['NMpri'], r['NMsupp'], r['MAPQpri'],
                                  r['MAPQsupp'], r['NP'], r['maxASsupp'], r['pe'] + r['supp'], r['pe'], r['supp'],
                                  r['sc'], r['block_edge'], r['raw_reads_10kb'], r['linked'], r['neigh'],
                                  r['ref_bases'], r['Prob']])
        else:
            format_fields[name] = [0] * 20

    return format_fields, cols



def to_vcf(df, args, names, outfile):

    outfile.write(HEADER + "\t" + "\t".join(names) + "\n")
    click.echo("Input samples: {}".format(str(list(names))), err=True)

    version = pkg_resources.require("dysgu")[0].version

    if len(names) > 1:
        dm = df.sort_values(["partners"], ascending=False)
    else:
        dm = df

    seen_idx = set([])

    count = 0
    recs = []
    for idx, r in dm.iterrows():

        if idx in seen_idx:
            continue

        format_f, df_rows = gen_format_fields(r, df, names)

        if "partners" in r:
            seen_idx |= set(r["partners"])

        r_main = make_main_record(r, version, count, format_f, df_rows, df)
        recs.append(r_main)
        count += 1


    for rec in sorted(recs, key=lambda x: (x[0], x[1])):

        outfile.write("\t".join(list(map(str, rec))) + "\n")


def to_csv(df, args, names, outfile):
    # Convert the partners to a human readable format
    keytable = ["chrA", "posA", "chrB", "posB", "sample", "id", "kind", "svtype", "join_type", "cipos95A", "cipos95B",
         "DP", "DN", "DApri", "DAsupp",  "NMpri", "NMsupp", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "pe", "supp", "sc", "block_edge",
         "raw_reads_10kb",
          "linked", "contigA", "contigB",  "gc", "neigh", "rep", "ref_bases", "Prob"]

    if "partners" not in df.columns:
        if "table_name" in df.columns:
            del df["table_name"]
        if "event_id" in df.columns:
            del df["event_id"]
        df[keytable].to_csv(outfile, index=False)
        return

    keytable.append("partners")

    if args["separate"] == "False":
        p2 = []
        for idx, r in df.iterrows():

            if len(r["partners"]) == 0:
                p2.append("")
                continue
            else:
                key = "|".join([f"{df.iloc[idx]['table_name']},{idx}" for idx in r["partners"]])
                p2.append(key)

        df["partners"] = p2

    elif args["separate"] == "True":

        for k, df2 in df.groupby("table_name"):
            df2 = df2.copy()
            # Convert parnet indexes into original indexes
            ori = []
            for item in df2["partners"]:
                if not item:
                    ori.append("")
                else:
                    ori.append("|".join([f"{df.iloc[i]['table_name']},{df.iloc[i]['id']}" for i in item]))
            df2["partners"] = ori

            del df2["table_name"]
            if "event_id" in df:
                del df["event_id"]
            df2[keytable].to_csv(outfile[k], index=False)

    else:
        del df["table_name"]
        del df["event_id"]
        df[keytable].to_csv(outfile, index=False)


def view_file(args):
    t0 = time.time()
    if args["separate"] == "True":
        if args["out_format"] == "vcf":
            raise ValueError("--separate only supported for --out-format csv")

        if not all(os.path.splitext(i)[1] == ".csv" for i in args["input_files"]):
            raise ValueError("All input files must have .csv extension")
    #
    # if args["out_format"]

    seen_names = sortedcontainers.SortedSet([])
    dfs = []
    # indexes = []
    for item in args["input_files"]:

        if item == "-":
            df = pd.read_csv(sys.stdin)
            name = list(set(df["sample"]))

            if len(name) > 1:
                raise ValueError("More than one sample in stdin")
            else:

                name = name[0]
                seen_names.add(name)
        else:
            name, ext = os.path.splitext(item)
            if ext != ".csv":
                raise ValueError("View accepts files with .csv extension")

            if name in seen_names:
                raise ValueError("Duplicate input file names")
            seen_names.add(name)
            df = pd.read_csv(item, index_col=None)

        if len(set(df["sample"])) > 1:
            raise ValueError("More than one sample per input file")

        df["contigA"] = [None if pd.isna(i) else i for i in df["contigA"]]
        df["contigB"] = [None if pd.isna(i) else i for i in df["contigB"]]

        df["table_name"] = [name for _ in range(len(df))]

        if args["merge_within"] == "True":
            df = merge_df(df, {}, merge_within_sample=True)

        if "partners" in df:
            del df["partners"]

        # indexes += list(df.index)
        dfs.append(df)

    df = pd.concat(dfs)
    # df["old_indexes"] = indexes
    if args["merge_across"] == "True":
        if len(seen_names) > 1:
            df = merge_df(df, seen_names)
            # echo("here", df[["chrA", "posA", "chrB", "posB", "partners"]])
    outfile = open_outfile(args)

    if args["out_format"] == "vcf":
        to_vcf(df, args, seen_names, outfile)
    else:
        to_csv(df, args, seen_names, outfile)

    click.echo("dysgu view complete h:m:s, {}".format(str(datetime.timedelta(seconds=int(time.time() - t0))),
                                                    time.time() - t0),
               err=True)
