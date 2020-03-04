#cython: language_level = 3

import os
import pandas as pd
import click
import sys
import sortedcontainers
import pkg_resources
import time
import datetime
from collections import defaultdict
from dysgu import io_funcs, cluster


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
        found = cluster.merge_events(potential, 25, tree, try_rev=False, pick_best=True, add_partners=True)
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
        found = cluster.merge_events(potential, 25, tree, try_rev=False, pick_best=True, add_partners=False)
        return pd.DataFrame.from_records(found)


def to_csv(df, args, names, outfile):
    # Convert the partners to a human readable format

    keytable = io_funcs.col_names()

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

    seen_names = sortedcontainers.SortedSet([])
    dfs = []

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

            df = pd.read_csv(item, index_col=None)
            name = list(set(df["sample"]))

            if len(name) > 1:
                raise ValueError("More than one sample in stdin")
            else:

                name = name[0]
                if name in seen_names:
                    raise IOError("Sample name is present in more than one input file {}".format(name))
                seen_names.add(name)


        if len(set(df["sample"])) > 1:
            raise ValueError("More than one sample per input file")

        if args["no_chr"] == "True":
            df["chrA"] = [i.replace("chr", "") for i in df["chrA"]]
            df["chrB"] = [i.replace("chr", "") for i in df["chrB"]]

        if args["no_contigs"] == "True":
            df["contigA"] = [None] * len(df)
            df["contigB"] = [None] * len(df)
        else:
            df["contigA"] = [None if pd.isna(i) else i for i in df["contigA"]]
            df["contigB"] = [None if pd.isna(i) else i for i in df["contigB"]]

        df["table_name"] = [name for _ in range(len(df))]

        if args["merge_within"] == "True":
            l_before = len(df)
            df = merge_df(df, {}, merge_within_sample=True)
            click.echo("{} rows before merge {}, rows after {}".format(name, l_before, len(df)), err=True)

        if "partners" in df:
            del df["partners"]

        # indexes += list(df.index)
        dfs.append(df)

    df = pd.concat(dfs)
    # df["old_indexes"] = indexes
    if args["merge_across"] == "True":
        if len(seen_names) > 1:
            df = merge_df(df, {}) #seen_names)
            # echo("here", df[["chrA", "posA", "chrB", "posB", "partners"]])

    df = df.sort_values(["chrA", "posA", "chrB", "posB"])
    outfile = open_outfile(args)

    if args["out_format"] == "vcf":
        io_funcs.to_vcf(df, args, seen_names, outfile)
    else:
        to_csv(df, args, seen_names, outfile)

    click.echo("dysgu view complete h:m:s, {}".format(str(datetime.timedelta(seconds=int(time.time() - t0))),
                                                    time.time() - t0),
               err=True)
