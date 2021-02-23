#cython: language_level = 3

import os
import pandas as pd
import sys
import sortedcontainers
import logging
import time
import datetime
from collections import defaultdict
from dysgu import io_funcs, cluster
from dysgu.map_set_utils import echo, EventResult


def open_outfile(args, names_dict):
    if args["separate"] == "True":
        outfiles = {}
        # for item in args["input_files"]:
        for item, name in names_dict.items():
            outname = f"{name}.{args['post_fix']}.csv"
            outfiles[name] = open(outname, "w")
        return outfiles
    if args["svs_out"] == "-" or args["svs_out"] is None:
        logging.info("SVs output to stdout")
        outfile = sys.stdout
    else:
        logging.info("SVs output to {}".format(args["svs_out"]))
        outfile = open(args["svs_out"], "w")
    return outfile


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    # https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def set_numeric(d):
    numeric = ['svlen_precise', 'n_expansion', 'stride', 'ref_poly_bases', 'ref_rep', 'compress',
               'su', 'pe', 'supp', 'sc', 'NP', 'maxASsupp', 'plus', 'minus', 'spanning', 'double_clips',
               'n_unmapped_mates', 'block_edge',
               'n_small_tlen', 'bnd', 'ras', 'fas', 'cipos95A', 'cipos95B']
    for k in numeric:
        if k not in d or d[k] is None:
            d[k] = 0
    return d


def merge_df(df, n_samples, merge_dist, tree=None, merge_within_sample=False):
    logging.info("Merge distance: {} bp".format(merge_dist))
    df.reset_index(inplace=True)
    df["event_id"] = df.index
    df["contig"] = df["contigA"]
    df["contig2"] = df["contigB"]
    # Assume:
    df["preciseA"] = [1] * len(df)
    df["preciseB"] = [1] * len(df)
    potential = [dotdict(set_numeric(i)) for i in df.to_dict("records")]

    # potential = [EventResult().from_dict(i) for i in potential]

    bad_i = set([])  # These could not be merged at sample level, SVs probably too close?
    if not merge_within_sample:
        found = cluster.merge_events(potential, merge_dist, tree, try_rev=False, pick_best=False, add_partners=True,
                                     aggressive_ins_merge=True,
                                     same_sample=False)
        ff = defaultdict(set)

        for f in found:
            if f.partners is None:
                ff[f.event_id] = set([])
                # ff[f["event_id"]] = []
            else:
                # Remove partners from same sample
                current = f["table_name"]
                targets = set([])
                passed = True
                for item in f["partners"]:  # Only merge with one row per sample
                    t_name = df.loc[item]["table_name"]
                    if t_name != current and t_name not in targets and len(targets) < n_samples:
                        #ff[f["event_id"]].add(item)
                        #ff[item].add(f["event_id"])
                        targets.add(t_name)
                    else:  # Merged with self event. Can happen with clusters of SVs with small spacing
                        bad_i.add(item)
                        passed = False
                if passed:  # enumerate support between components
                    g = f["partners"] + [f["event_id"]]
                    for t1 in g:
                        for t2 in g:
                            if t2 != t1:
                                ff[t1].add(t2)

        df = df.drop(bad_i)
        df["partners"] = [ff[i] if i in ff else set([]) for i in df.index]
        return df
    else:
        found = cluster.merge_events(potential, merge_dist, tree, try_rev=False, pick_best=True, add_partners=False,
                                     same_sample=True)
        return pd.DataFrame.from_records(found)


def to_csv(df, args, names, outfile, extended, small_output):
    # Convert the partners to a human readable format
    keytable = io_funcs.col_names(extended, small_output)
    fmt = keytable.pop()
    keytable += fmt
    # keytable = [item for sublist in keytable for item in sublist]
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
                key = "|".join([f"{df.loc[idx]['table_name']},{idx}" for idx in r["partners"]])
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
                    ori.append("|".join([f"{df.loc[i]['table_name']},{df.loc[i]['id']}" for i in item]))
            df2["partners"] = ori
            del df2["table_name"]
            if "event_id" in df:
                del df["event_id"]
            df2[keytable].to_csv(outfile[k], index=False)
    else:
        del df["table_name"]
        del df["event_id"]
        df[keytable].to_csv(outfile, index=False)


def vcf_to_df(path):

    # Parse sample
    header = ""
    with open(path, "r") as fin:
        last = ""
        for line in fin:
            if line[0] == "#":
                if last:
                    header += last
                last = line
            else:
                header += "\t".join(last.split("\t")[:9])
                break
        sample = last.strip().split("\t")[9]

    df = pd.read_csv(path, index_col=None, comment="#", sep="\t", header=None)
    if len(df.columns) > 10:
        raise ValueError(f"Can only merge files with one sample in. N samples = {len(df.columns) - 9}")
    parsed = pd.DataFrame()
    parsed["chrA"] = df[0]
    parsed["posA"] = df[1]
    parsed["id"] = df[2]
    parsed["sample"] = [sample] * len(df)
    info = []
    for k in list(df[7]):
        if k:
            info.append(dict(i.split("=") for i in k.split(";") if "=" in i))

    n_fields = None
    for idx, (k1, k2), in enumerate(zip(df[8], df[9])):  # Overwrite info column with anything in format
        if n_fields is None:
            n_fields = len(k1.split(":"))
        info[idx].update({i: j for i, j in zip(k1.split(":"), k2.split(":"))})

    info_df = pd.DataFrame.from_records(info)
    df = pd.concat([parsed, info_df], axis=1)

    col_map = {"chrA": ("chrA", str),
               "posA": ("posA", int),
               "CHR2": ("chrB", str),
               "GRP": ("grp_id", int),
               "NGRP": ("n_in_grp", int),
               "END": ("posB", int),
               "sample": ("sample", str),
               "id": ("id", str),
               "KIND": ("kind", str),
               "SVTYPE": ("svtype", str),
               "CT": ("join_type", str),
               "CIPOS95": ("cipos95A", int),
               "CIEND95": ("cipos95B", int),
               "NMP": ("NMpri", float),
               "NMB": ("NMbase", float),
               "NMS": ("NMsupp", float),
               "MAPQP": ("MAPQpri", float),
               "MAPQS": ("MAPQsupp", float),
               "NP": ("NP", int),
               "OL": ("query_overlap", int),
               "MAS": ("maxASsupp", int),
               "SU": ("su", int),
               "WR": ("spanning", int),
               "PE": ("pe", int),
               "SR": ("supp", int),
               "SC": ("sc", int),
               "BND": ("bnd", int),
               "SQC": ("sqc", float),
               "SCW": ("scw", float),
               "RT": ("type", str),
               "BE": ("block_edge", int),
               "COV": ("raw_reads_10kb", float),
               "MCOV": ("mcov", float),
               "LNK": ("linked", int),
               "CONTIGA": ("contigA", str),
               "CONTIGB": ("contigB", str), "GC": ("gc", float),
               "NEIGH": ("neigh", int),
               "NEIGH10": ("neigh10kb", int),
               "REP": ("rep", float),
               "REPSC": ("rep_sc", float),
               "LPREC": ("svlen_precise", int),
               "NEXP": ("n_expansion", int),
               "STRIDE": ("stride", int),
               "EXPSEQ": ("exp_seq", str),
               "RPOLY": ("ref_poly_bases", int),
               "GT": ("GT", str),
               "GQ": ("GQ", object),
               "RB": ("ref_bases", int),
               "SVLEN": ("svlen", int),
               "PS": ("plus", int),
               "MS": ("minus", int),
               "PROB": ("prob", float),
               "NG": ("n_gaps", float),
               "NSA": ("n_sa", float),
               "NXA": ("n_xa", float),
               "NMU": ("n_unmapped_mates", int),
               "NDC": ("double_clips", int),
               "RMS": ("remap_score", int),
               "RED": ("remap_ed", int),
               "BCC": ("bad_clip_count", int),
               "STL": ("n_small_tlen", int),
               "RAS": ("ras", int),
               "FAS": ("fas", int),
               "ICN": ("inner_cn", float),
               "OCN": ("outer_cn", float),
               "CMP": ("compress", float),
               "FCC": ("fcc", float),
               "RR": ("rep_ref", float),
               "JIT": ("jitter", float)
               }

    df.rename(columns={k: v[0] for k, v in col_map.items()}, inplace=True)

    for k, dtype in col_map.values():
        if k in df:
            if df[k].dtype != dtype:
                if dtype == str:
                    df[k] = df[k].fillna("")
                else:
                    df[k] = df[k].fillna(0)
                try:
                    df[k] = df[k].astype(dtype)
                except ValueError:
                    logging.info("Problem for feature {}, could not intepret as {}".format(k, dtype))
                    echo(list(df[k]))
                    quit()
    if "contigA" not in df:
        df["contigA"] = [""] * len(df)
    if "contigB" not in df:
        df["contigB"] = [""] * len(df)

    df["posA"] = df["posA"].astype(int)
    df["posB"] = df["posB"].astype(int)

    return df, header, n_fields


def view_file(args):

    t0 = time.time()
    if args["separate"] == "True":
        if args["out_format"] == "vcf":
            raise ValueError("--separate only supported for --out-format csv")

        if not all(os.path.splitext(i)[1] == ".csv" for i in args["input_files"]):
            raise ValueError("All input files must have .csv extension")

    args["metrics"] = False  # only option supported so far
    args["contigs"] = False
    seen_names = sortedcontainers.SortedSet([])
    names_dict = {}
    name_c = defaultdict(int)
    dfs = []

    header = None

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
            if ext != ".csv":  # assume vcf
                df, header, _ = vcf_to_df(item)  # header here, assume all input has same number of fields
            else:
                df = pd.read_csv(item, index_col=None)
            name = list(set(df["sample"]))

            if len(name) > 1:
                raise ValueError("More than one sample in stdin")
            else:
                name = name[0]
                names_dict[item] = name
                if name in seen_names:
                    logging.info("Sample {} is present in more than one input file".format(name))
                    bname = f"{name}_{name_c[name]}"
                    df["sample"] = [bname] * len(df)
                    df["table_name"] = [bname for _ in range(len(df))]
                    seen_names.add(bname)
                else:
                    df["table_name"] = [name for _ in range(len(df))]
                    seen_names.add(name)
                name_c[name] += 1

        if len(set(df["sample"])) > 1:
            raise ValueError("More than one sample per input file")

        if args["no_chr"] == "True":
            df["chrA"] = [i.replace("chr", "") for i in df["chrA"]]
            df["chrB"] = [i.replace("chr", "") for i in df["chrB"]]

        if args["no_contigs"] == "True":
            df["contigA"] = [None] * len(df)
            df["contigB"] = [None] * len(df)
        else:
            # click.echo(df.columns, err=True)
            df["contigA"] = [None if pd.isna(i) else i for i in df["contigA"]]
            df["contigB"] = [None if pd.isna(i) else i for i in df["contigB"]]

        if args["merge_within"] == "True":
            l_before = len(df)
            df = merge_df(df, 1, args["merge_dist"], {}, merge_within_sample=True)
            logging.info("{} rows before merge-within {}, rows after {}".format(name, l_before, len(df)))

        if "partners" in df:
            del df["partners"]
        dfs.append(df)

    df = pd.concat(dfs)
    if args["merge_across"] == "True":
        if len(seen_names) > 1:
            df = merge_df(df, len(seen_names), args["merge_dist"], {})

    df = df.sort_values(["chrA", "posA", "chrB", "posB"])
    outfile = open_outfile(args, names_dict)

    if args["out_format"] == "vcf":
        count = io_funcs.to_vcf(df, args, seen_names, outfile, header=header, extended_tags=False, small_output_f=True)
        logging.info("Sample rows before merge {}, rows after {}".format(list(map(len, dfs)), count))
    else:
        to_csv(df, args, seen_names, outfile, True, False)

    logging.info("dysgu merge complete h:m:s, {}".format(str(datetime.timedelta(seconds=int(time.time() - t0))),
                                                    time.time() - t0))
