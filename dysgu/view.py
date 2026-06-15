import glob
import os
import sys
import pandas as pd
import sortedcontainers
from superintervals import IntervalMap
import logging
import time
import datetime
from collections import defaultdict
from dysgu import io_funcs, cluster, merge_svs
from dysgu.map_set_utils import echo, merge_intervals
import multiprocessing
import queue
import gzip
import pysam
import numpy as np
from scipy.spatial.distance import cosine
import networkx as nx
from sys import stdout
from heapq import heappop, heappush
import resource
import re
import traceback


def open_outfile(args, names_list, log_messages=True):
    if args["separate"] == "True":
        outfiles = {}
        for name in names_list:
            outname = f"{name}.{args['post_fix']}.csv"
            outfiles[name] = open(outname, "w")
        return outfiles
    if args["svs_out"] == "-" or args["svs_out"] is None:
        if log_messages:
            logging.info("SVs output to stdout")
        outfile = stdout
    else:
        if log_messages:
            logging.info("SVs output to {}".format(args["svs_out"]))
        outfile = open(args["svs_out"], "w")
    return outfile


def _run_process_job(func, func_args, result_queue, job_index, job_label):
    try:
        result = func(*func_args)
    except BaseException:
        result_queue.put((job_index, job_label, "error", traceback.format_exc()))
    else:
        result_queue.put((job_index, job_label, "ok", result))


def run_process_jobs(func, job_args, procs, job_labels=None, progress=False):
    if not job_args:
        return []

    procs = max(1, min(procs, len(job_args)))
    job_labels = job_labels or [str(i) for i in range(len(job_args))]
    results = [None] * len(job_args)
    completed = set()
    failures = []
    result_queue = multiprocessing.Queue()
    running = {}
    next_job = 0
    n_jobs = len(job_args)
    finished = 0

    def start_job(job_index):
        label = job_labels[job_index]
        if progress:
            logging.info(f"[{job_index + 1}/{n_jobs}] started {label}")
        proc = multiprocessing.Process(
            target=_run_process_job,
            args=(func, job_args[job_index], result_queue, job_index, label),
        )
        proc.start()
        running[job_index] = proc

    def drain_results(block=False):
        nonlocal finished
        got_any = False
        while True:
            try:
                if block:
                    job_index, label, status, payload = result_queue.get(timeout=0.2)
                else:
                    job_index, label, status, payload = result_queue.get_nowait()
            except queue.Empty:
                break
            got_any = True
            completed.add(job_index)
            if status == "ok":
                results[job_index] = payload
                if progress:
                    finished += 1
                    logging.info(f"[{finished}/{n_jobs}] finished {label}")
            else:
                failures.append((label, payload))
            block = False
        return got_any

    while next_job < len(job_args) or running:
        while next_job < len(job_args) and len(running) < procs:
            start_job(next_job)
            next_job += 1

        got_result = drain_results()

        for job_index, proc in list(running.items()):
            if proc.exitcode is None:
                continue
            proc.join()
            del running[job_index]
            if job_index not in completed and proc.exitcode != 0:
                failures.append((job_labels[job_index], f"worker exited with code {proc.exitcode}"))

        if failures:
            for proc in running.values():
                proc.terminate()
            for proc in running.values():
                proc.join()
            msg = "\n\n".join(f"{label}:\n{err}" for label, err in failures)
            raise RuntimeError(f"One or more merge worker jobs failed:\n{msg}")

        if not got_result:
            time.sleep(0.1)

    while len(completed) < len(job_args) and drain_results(block=True):
        pass
    drain_results()

    missing = [job_labels[i] for i in range(len(job_args)) if i not in completed]
    if failures or missing:
        messages = [f"{label}:\n{err}" for label, err in failures]
        messages += [f"{label}:\nworker exited without returning a result" for label in missing]
        raise RuntimeError("One or more merge worker jobs failed:\n" + "\n\n".join(messages))

    return results


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


def mung_df(df, args):
    if args["verbosity"] == "0":
        df["contigA"] = [None] * len(df)
        df["contigB"] = [None] * len(df)
    else:
        df["contigA"] = [None if pd.isna(i) else i for i in df["contigA"]]
        df["contigB"] = [None if pd.isna(i) else i for i in df["contigB"]]
    return df


def merge_across_samples(df, potential, merge_dist, tree, aggressive, samples, progressive):

    samples = list(samples)
    # logging.info(f"{samples[0]}")

    if progressive:
        by_sample = defaultdict(list)
        for e in potential:
            by_sample[e.sample].append(e)

        d1 = by_sample.get(samples[0], [])
        components = {e.event_id: {e.event_id} for e in d1}
        for idx, samp in enumerate(samples[1:], start=1):
            d2 = by_sample.get(samp, [])
            pot = d1 + d2

            found = merge_svs.merge_events(pot, merge_dist, tree, try_rev=False, pick_best=False, add_partners=True,
                                           aggressive_ins_merge=True, same_sample=False, procs=1)

            new_components = {}
            for f in found:
                if f.partners is None:
                    new_components[f.event_id] = components.get(f.event_id, {f.event_id})
                else:
                    members = set()
                    for item in f["partners"] + [f["event_id"]]:
                        members.update(components.get(item, {item}))
                    new_components[f.event_id] = members

            d1 = found
            components = new_components
            # logging.info(f"Merged {samp} ({idx}/{len(samples)}), cohort SV sites: {len(found)}")

    else:
        found = merge_svs.merge_events(potential, merge_dist, tree, try_rev=False, pick_best=False, add_partners=True,
                                       aggressive_ins_merge=True, same_sample=False)
        components = {}

        for f in found:
            if f.partners is None:
                continue
            else:
                # Remove partners from same sample
                current = f["table_name"]
                targets = set([])
                passed = True
                for item in f["partners"]:  # Only merge with one row per sample
                    t_name = df.loc[item]["table_name"]
                    if t_name != current and t_name not in targets and len(targets) < len(samples):
                        targets.add(t_name)

                    elif not aggressive:
                        # Merged with self event. Can happen with clusters of SVs with small spacing
                        # e.g. a is merged with b and c, where a is from sample1 and b and c are from sample2
                        # safer not to merge? otherwise variants can be lost
                        passed = False
                if passed:  # enumerate support between components
                    components[f.event_id] = set(f["partners"] + [f["event_id"]])

    partner_map = {}
    for rep, members in components.items():
        if len(members) > 1:
            members.discard(rep)
            partner_map[rep] = members
    skip_partner = pd.Series(False, index=df.index)
    for members in partner_map.values():
        skip_partner.loc[list(members)] = True
    df["partners"] = [partner_map.get(i) for i in df.index]
    df["_skip_partner"] = skip_partner
    return df


def merge_df(df, samples, merge_dist, tree=None, merge_within_sample=False, aggressive=False, log_messages=True,
             progressive=False):
    if log_messages:
        logging.info("Merge distance: {} bp".format(merge_dist))
    df.reset_index(inplace=True)
    df["event_id"] = df.index
    df["contig"] = df["contigA"]
    df["contig2"] = df["contigB"]
    # Assume:
    df["preciseA"] = [1] * len(df)
    df["preciseB"] = [1] * len(df)
    potential = [dotdict(set_numeric(i)) for i in df.to_dict("records")]

    if not merge_within_sample:
        return merge_across_samples(df, potential, merge_dist, tree, aggressive, samples, progressive)
    else:
        found = merge_svs.merge_events(potential, merge_dist, tree, try_rev=False, pick_best=True, add_partners=False,
                                     same_sample=True, aggressive_ins_merge=True)
        return pd.DataFrame.from_records(found)


def to_csv(df, args, outfile, small_output):
    keytable = io_funcs.col_names(small_output)
    fmt = keytable.pop()
    keytable += fmt
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
            partners = r["partners"]
            partner_collection = isinstance(partners, (set, list, tuple))
            if partners is None or (partner_collection and len(partners) == 0) or (
                    not partner_collection and (pd.isna(partners) or len(partners) == 0)):
                p2.append("")
                continue
            else:
                key = "|".join([f"{df.at[idx, 'table_name']},{idx}" for idx in r["partners"]])
                p2.append(key)
        df["partners"] = p2
        df[keytable].to_csv(outfile, index=False)

    elif args["separate"] == "True":
        for k, df2 in df.groupby("table_name"):
            df2 = df2.copy()
            ori = []
            for item in df2["partners"]:
                if not item:
                    ori.append("")
                else:
                    ori.append("|".join([f"{df.at[i, 'table_name']},{df.at[i, 'event_id']}" for i in item]))
            df2["partners"] = ori
            del df2["table_name"]
            df2[keytable].to_csv(outfile[k], index=False)
    else:
        del df["table_name"]
        del df["event_id"]
        df[keytable].to_csv(outfile, index=False)


def read_from_inputfile(path):
    if path != '-' and isinstance(path, str):
        try:
            with open(path, "r") as fin:
                for line in fin:
                    yield line
        except UnicodeDecodeError:
            with gzip.open(path, 'r') as fin:
                for line in fin:
                    yield line.decode('ascii')

    elif path == '-':
        for line in sys.stdin:
            yield line
    else:
        for line in path:
            yield line


def _parse_vcf_columns(df, samples):
    """Parse a DataFrame of raw VCF columns (0-8 = standard fields, 9 = the single FORMAT-value
    column for that record) into dysgu's internal feature DataFrame.

    `samples` is a per-row sample-name sequence (length == len(df)). For a single-sample input every
    row shares the same name; for an internal blob file each row carries its own originating sample.
    Returns (df, n_fields).
    """
    parsed = pd.DataFrame()
    parsed["chrA"] = df[0]
    parsed["posA"] = df[1]
    parsed["event_id"] = df[2]
    parsed["ref_seq"] = df[3]
    parsed["variant_seq"] = df[4]
    parsed["filter"] = df[6]
    parsed["sample"] = list(samples)
    info = []
    for k in list(df[7]):
        if k:
            info.append(dict(i.split("=") for i in k.split(";") if "=" in i))
        else:
            info.append({})

    fmt_cols = df[8].str.split(':')
    sample_cols = df[9].str.split(':')
    n_fields = None
    for idx, (k1, k2) in enumerate(zip(fmt_cols, sample_cols)):  # Overwrite info column with anything in format
        if n_fields is None:
            n_fields = len(k1)
        info[idx].update({i: j for i, j in zip(k1, k2)})

    info_df = pd.DataFrame.from_records(info)
    df = pd.concat([parsed, info_df], axis=1)

    col_map = {"chrA": ("chrA", str),
               "posA": ("posA", np.int64),
               "CHR2": ("chrB", str),
               "GRP": ("grp_id", np.int64),
               "NGRP": ("n_in_grp", np.int64),
               "END": ("posB", np.int64),
               "CHR2_POS": ("posB_tra", np.int64),
               "sample": ("sample", str),
               "event_id": ("event_id", np.int64),
               "KIND": ("kind", str),
               "SVTYPE": ("svtype", str),
               "CT": ("join_type", str),
               "CIPOS95": ("cipos95A", np.int64),
               "CIEND95": ("cipos95B", np.int64),
               "NMP": ("NMpri", np.float64),
               "NMB": ("NMbase", np.float64),
               "NMS": ("NMsupp", np.float64),
               "MAPQP": ("MAPQpri", np.float64),
               "MAPQS": ("MAPQsupp", np.float64),
               "NP": ("NP", np.int64),
               "OL": ("query_overlap", np.int64),
               "MAS": ("maxASsupp", np.int64),
               "SU": ("su", np.int64),
               "WR": ("spanning", np.int64),
               "PE": ("pe", np.int64),
               "SR": ("supp", np.int64),
               "SC": ("sc", np.int64),
               "BND": ("bnd", np.int64),
               "SQC": ("sqc", np.float64),
               "SCW": ("scw", np.float64),
               "SQR": ("clip_qual_ratio", np.float64),
               "RT": ("type", str),
               "BE": ("block_edge", np.int64),
               "COV": ("raw_reads_10kb", np.float64),
               "MCOV": ("mcov", np.float64),
               "LNK": ("linked", np.int64),
               "CONTIGA": ("contigA", str),
               "CONTIGB": ("contigB", str),
               "ref_seq": ("ref_seq", str),
               "variant_seq": ("variant_seq", str),
               "GC": ("gc", np.float64),
               "NEIGH": ("neigh", np.int64),
               "NEIGH10": ("neigh10kb", np.int64),
               "REP": ("rep", np.float64),
               "REPSC": ("rep_sc", np.float64),
               "LPREC": ("svlen_precise", np.int64),
               "NEXP": ("n_expansion", np.int64),
               "STRIDE": ("stride", np.int64),
               "EXPSEQ": ("exp_seq", str),
               "RPOLY": ("ref_poly_bases", np.int64),
               "GT": ("GT", str),
               "GQ": ("GQ", object),
               "RB": ("ref_bases", np.int64),
               "SVLEN": ("svlen", np.int64),
               "PS": ("plus", np.int64),
               "MS": ("minus", np.int64),
               "SBT": ("strand_binom_t", np.float64),
               "PROB": ("prob", np.float64),
               "NG": ("n_gaps", np.float64),
               "NSA": ("n_sa", np.float64),
               "NXA": ("n_xa", np.float64),
               "NMU": ("n_unmapped_mates", np.int64),
               "NDC": ("double_clips", np.int64),
               "RMS": ("remap_score", np.int64),
               "RED": ("remap_ed", np.int64),
               "BCC": ("bad_clip_count", np.int64),
               "STL": ("n_small_tlen", np.int64),
               "RAS": ("ras", np.int64),
               "FAS": ("fas", np.int64),
               "ICN": ("inner_cn", np.float64),
               "OCN": ("outer_cn", np.float64),
               "CMP": ("compress", np.float64),
               "FCC": ("fcc", np.float64),
               "RR": ("ref_rep", np.float64),
               "JIT": ("jitter", np.float64),
               "LEFT_SVINSSEQ": ("left_ins_seq", str),
               "RIGHT_SVINSSEQ": ("right_ins_seq", str),
               "PSET": ("phase_set", np.int64),
               "HP": ("haplotype", str),
               "AF": ("a_freq", np.float64),
               }

    # First check which original columns are missing before renaming
    original_columns_in_df = set(df.columns)
    missing_original_columns = {}

    # Required columns with default values:
    required = {"phase_set": -1, "haplotype": "-1", "a_freq": -1, "posB_tra": -1, "svlen": -1, "MAPQsupp": -1,
                "NMbase": -1, "exp_seq": ""}
    for orig_col, (new_col, _) in col_map.items():
        if orig_col not in original_columns_in_df and new_col in required:
            missing_original_columns[new_col] = required[new_col]

    # Now rename the columns that are present
    df.rename(columns={k: v[0] for k, v in col_map.items() if k in df.columns}, inplace=True)

    # Add the missing columns after renaming
    for value_name, default_value in missing_original_columns.items():
        df[value_name] = [default_value] * len(df)

    df["GQ"] = pd.to_numeric(df["GQ"], errors='coerce').fillna(".")
    for k, dtype in col_map.values():
        if k in df:
            if df[k].dtype != dtype:
                if dtype == str:
                    if k == "GT":
                        df[k] = ["0/0" if not i else i for i in list(df[k])]
                    else:
                        df[k] = df[k].fillna("")
                else:
                    # VCF missing value "." can appear in numeric INFO/FORMAT fields
                    # (e.g. for 0/0 records emitted by --all-sites True); coerce to NaN then fill.
                    df[k] = pd.to_numeric(df[k], errors='coerce').fillna(0)
                try:
                    df[k] = df[k].astype(dtype)
                except (ValueError, OverflowError) as e:
                    raise ValueError("Problem for feature {}, could not interpret as {}".format(k, dtype)) from e
    if "contigA" not in df:
        df["contigA"] = [""] * len(df)
    if "contigB" not in df:
        df["contigB"] = [""] * len(df)
    if "posB" not in df:
        if "posB_tra" in df:
            df["posB"] = df["posB_tra"]
        else:
            df["posB"] = df["posA"]
    if 'posB_tra' in df:
        df["posB"] = [i if svt != 'TRA' else j for i, j, svt in zip(df['posB'], df['posB_tra'], df['svtype'])]
        del df['posB_tra']
    df["posA"] = df["posA"].astype(int) - 1  # convert to 0-indexing
    df["posB"] = df["posB"].astype(int) - 1
    return df, n_fields


def _read_vcf_header(path):
    """Read the ## header block of a VCF, returning (header_str, contig_names_str, first_data_line).

    header_str ends with the truncated #CHROM..FORMAT (first 9) columns of the data section, matching
    the original vcf_to_df behaviour.
    """
    fin = read_from_inputfile(path)
    header = ""
    last = ""
    contig_names = ""
    for line in fin:
        if line[:2] == "##":
            header += line
            if line.startswith("##contig="):
                contig_names += line
        else:
            header += "\t".join(line.split("\t")[:9])
            last = line
            break
    return header, contig_names, last


def vcf_to_df(path):
    if os.stat(path).st_size == 0:
        logging.warning(f"File is empty: {path}")
        return pd.DataFrame()

    header, contig_names, last = _read_vcf_header(path)
    if not header:
        logging.critical(f"File does not have vcf header: {path}")
        quit()

    sample = last.strip().split("\t")[9]
    if path == '-':
        path_h = sys.stdin
    else:
        path_h = path
    df = pd.read_csv(path_h, index_col=None, comment="#", sep="\t", header=None, dtype=str)
    if len(df.columns) > 10:
        msg = f"Can only merge files with one sample in. N samples = {len(df.columns) - 9}"
        raise ValueError(msg)
    df, n_fields = _parse_vcf_columns(df, [sample] * len(df))
    return df, header, n_fields, "\n" + contig_names.strip() if contig_names else contig_names


def get_names_list(file_list, ignore_csv=True):
    seen_names = sortedcontainers.SortedSet([])
    names_list = []
    name_c = defaultdict(int)
    for item in file_list:
        if item == "-":
            raise ValueError("Reading from stdin is not supported using merge")
        if item.endswith(".csv"):
            if ignore_csv:
                continue
            raise ValueError(".csv files are not supported when using merge and option --wd")
        name = list(pysam.VariantFile(item, 'r').header.samples)
        if len(name) > 1:
            msg = f"Sample file {item} contains more than one sample: {name}"
            raise ValueError(msg)
        name = name[0]
        if name in seen_names:
            bname = f"{name}_{name_c[name]}"
            names_list.append(bname)
            logging.info(f"Sample {name} is present in more than one input file, sample in file {item} will be named {bname} in merged vcf")
            seen_names.add(bname)
        else:
            names_list.append(name)
            seen_names.add(name)
        name_c[name] += 1
    return seen_names, names_list


def blob_to_df(path):
    """Read an internal blob file into dysgu's feature DataFrame.

    A blob line is a tab-separated record: the originating sample name followed by the 10 columns of
    that sample's dysgu VCF record (CHROM..FORMAT, then the single sample's FORMAT values). Records
    from many samples are mixed in one file; the per-row sample column preserves their identity so
    merge-across works exactly as it did with one-file-per-sample input.
    """
    if os.stat(path).st_size == 0:
        return pd.DataFrame()
    raw = pd.read_csv(path, index_col=None, sep="\t", header=None, dtype=str)
    samples = raw[0]
    # Shift columns left by one so column i holds standard VCF field i (0..9), matching vcf_to_df.
    body = raw.iloc[:, 1:]
    body.columns = range(body.shape[1])
    df, _ = _parse_vcf_columns(body, samples)
    return df


def process_file_list(args, file_list, seen_names, names_list, log_messages, show_progress=False, job_id=None,
                      header=None, contig_names=None):
    """Merge one job's records and write the intermediate VCF.

    file_list is a list of blob files (each may hold records from many samples, one sample per row).
    header/contig_names provide the VCF header for output (captured once from an input file) since a
    blob carries no usable single header of its own. The legacy one-file-per-sample path also routes
    here, in which case each file is a single-sample VCF/CSV read via vcf_to_df.
    """
    dfs = []
    for fi, item in enumerate(file_list):
        _, ext = os.path.splitext(item)
        if ext == ".csv":
            df = pd.read_csv(item, index_col=None)
        elif ext == ".blob":
            df = blob_to_df(item)
        else:  # single-sample vcf (legacy path with no working directory): one bname per file
            df, header, _, contig_names = vcf_to_df(item)
            df["sample"] = [names_list[fi]] * len(df)
        if not len(df):
            continue
        df["table_name"] = list(df["sample"])
        df = mung_df(df, args)
        dfs.append(df)

    if not dfs:
        # Still emit an (empty-bodied) header-only file so downstream assembly has a valid VCF.
        outfile = open_outfile(args, names_list, log_messages=log_messages)
        try:
            if args["out_format"] == "vcf" and header is not None:
                io_funcs.to_vcf(pd.DataFrame(), args, seen_names, outfile, header=header, small_output_f=True,
                                contig_names=contig_names, show_names=log_messages)
        finally:
            if outfile is not stdout:
                outfile.close()
        return

    df = pd.concat(dfs)

    if args["merge_within"] == "True":
        # merge-within is per-sample: run it on each sample's records independently.
        l_before = len(df)
        if show_progress and job_id:
            logging.info("{} started, records {}".format(job_id, l_before))
        merged_parts = []
        for _, sub in df.groupby("table_name", sort=False):
            sub = merge_df(sub.copy(), 1, args["merge_dist"], {}, merge_within_sample=True,
                           aggressive=args['collapse_nearby'] == "True", log_messages=False,
                           progressive=args['progressive'])
            if "partners" in sub:
                del sub["partners"]
            merged_parts.append(sub)
        df = pd.concat(merged_parts)
        if show_progress and job_id:
            logging.info("{} rows before merge-within {}, rows after {}".format(job_id, l_before, len(df)))

    if args["merge_across"] == "True" and len(seen_names) > 1:
        df = merge_df(df, seen_names, args["merge_dist"], {}, aggressive=args['collapse_nearby'] == "True",
                      log_messages=log_messages, progressive=args['progressive'])

    df = df.sort_values(["chrA", "posA", "chrB", "posB"])

    outfile = open_outfile(args, names_list, log_messages=log_messages)
    try:
        if args["out_format"] == "vcf":
            if show_progress and job_id:
                logging.info("{} started, records {}".format(job_id, len(df)))
            count = io_funcs.to_vcf(df, args, seen_names, outfile, header=header, small_output_f=True,
                                    contig_names=contig_names, show_names=log_messages)
            if log_messages:
                logging.info("Rows after merge {}".format(count))
            elif show_progress and job_id:
                logging.info("{} rows after {}".format(job_id, count))
        else:
            to_csv(df, args, outfile, small_output=False)
    finally:
        if outfile is not stdout:
            if isinstance(outfile, dict):
                for out in outfile.values():
                    out.close()
            else:
                outfile.close()


class VcfWriter(object):
    def __init__(self, out_path, target_header, new_name=None):
        self.path = out_path
        self.vcf = open(out_path, 'w') if out_path != '-' else stdout
        self.header = target_header
        str_hdr = str(target_header)
        if new_name:
            endi = len(str_hdr)
            while True:
                endi -= 1
                if str_hdr[endi] == "\t":
                    break
            str_hdr = str_hdr[:endi + 1] + f"{new_name}\n"
        self.vcf.write(str_hdr)
    def close(self):
        self.vcf.close()
    def write(self, r):
        self.vcf.write(str(r))


def _record_group_key(r):
    """Return (group_key, anchor_pos): the grouping two SVs must share to possibly merge, plus the
    canonical anchor position used to bin the record into a blob.

    INS and DUP collapse to INSDUP. For intra-chromosomal SVs the group is the chromosome and the
    anchor is the record's own position. For TRA (inter-chromosomal) the group is the sorted
    chromosome pair; merge_events canonicalises these so the merge anchor is the breakpoint on the
    lexicographically *smaller* chromosome. We therefore bin a TRA by that canonical position, so the
    two possible orientations of the same translocation (emitted as chrA->chrB in one sample and
    chrB->chrA in another) bin to the same coordinate and land in the same blob.
    """
    svtype = r.info['SVTYPE']
    svtype = 'INSDUP' if (svtype == "INS" or svtype == "DUP") else svtype
    if svtype != "TRA":
        return (svtype, r.chrom), r.pos
    chrom2 = r.info["CHR2"]
    anchor = r.pos if r.chrom <= chrom2 else int(r.info["CHR2_POS"])
    return (svtype, f'{min(r.chrom, chrom2)}_{max(r.chrom, chrom2)}'), anchor


def scan_sample_keys(item_path, name):
    """Pass 1 worker: scan one input VCF and return its grouping keys.

    Returns (name, keys, n_rows) where keys is a list of (group_key, anchor_pos). Order does not
    matter - boundaries are computed from globally sorted anchors, and each record is binned to its
    blob independently, so no pre-sorted input is required.
    """
    vcf = pysam.VariantFile(item_path, 'r')
    keys = []
    rows = 0
    for r in vcf.fetch():
        keys.append(_record_group_key(r))
        rows += 1
    vcf.close()
    return name, keys, rows


def compute_blob_boundaries(per_sample_keys, merge_dist, target_rows):
    """Decide blob boundaries from every sample's grouping keys.

    All records are placed in one global order by (group_index, anchor_pos). A blob is a contiguous
    run of ~target_rows records in that order. To keep blobs full even when individual chromosomes or
    SV types are small, a blob is allowed to span *multiple* groups - the merge step handles mixed
    groups in one file, and records in different groups can never merge anyway.

    A blob is only ended at a *safe* boundary so no mergeable pair is split: either a group change
    (records on different chromosomes / SV types never merge) or a within-group gap wider than
    2 * merge_dist. We don't *force* a cut at every group change - only once the blob has reached
    target_rows.

    Returns:
      group_order: dict group_key -> its global ordering index (groups in first-seen order).
      boundaries:  sorted list of (group_index, start_pos) blob-start keys. A record maps to the last
                   blob whose start key is <= its own (group_index, anchor_pos).
    """
    group_order = {}
    group_positions = defaultdict(list)
    for keys in per_sample_keys:
        for gk, pos in keys:
            if gk not in group_order:
                group_order[gk] = len(group_order)
            group_positions[gk].append(pos)

    gap_threshold = 2 * int(merge_dist)
    ordered_groups = sorted(group_order, key=lambda g: group_order[g])

    # Flatten into one globally ordered stream of (group_index, pos).
    stream = []
    for gk in ordered_groups:
        gi = group_order[gk]
        positions = group_positions[gk]
        positions.sort()
        for p in positions:
            stream.append((gi, p))

    n = len(stream)
    boundaries = []
    if n == 0:
        return group_order, boundaries

    boundaries.append(stream[0])  # first blob starts at the first record
    seg_count = 0
    for i in range(n):
        seg_count += 1
        if i + 1 < n:
            cur_gi, cur_pos = stream[i]
            nxt_gi, nxt_pos = stream[i + 1]
            # A safe place to start a new blob: different group, or a large within-group gap.
            safe = (nxt_gi != cur_gi) or (nxt_pos - cur_pos > gap_threshold)
            if seg_count >= target_rows and safe:
                boundaries.append(stream[i + 1])
                seg_count = 0
    return group_order, boundaries


def _blob_index_for(group_order, boundaries, gk, pos):
    """Map a record's (group_key, pos) to its blob index via binary search over blob-start keys.

    boundaries is the sorted list of (group_index, start_pos) starts; the owning blob is the last one
    whose start key is <= (this record's group_index, pos).
    """
    target = (group_order[gk], pos)
    lo, hi = 0, len(boundaries)
    while lo < hi:
        mid = (lo + hi) // 2
        if boundaries[mid] <= target:
            lo = mid + 1
        else:
            hi = mid
    return lo - 1


def write_blobs(wd, input_files, names_list, group_order, boundaries, show_progress):
    """Pass 2 (serial): stream every input VCF and write each record to its blob file.

    Each blob line is: <sample_name>\\t<original 10-column VCF record>. Running single-threaded keeps
    one writer per blob with no locking; writing is I/O-bound while the expensive merge stays
    parallel. Returns (blob_paths_in_order, input_rows).
    """
    blob_paths = [os.path.join(wd, f"blob_{i}.blob") for i in range(len(boundaries))]
    handles = [open(p, 'w') for p in blob_paths]
    input_rows = {}
    try:
        for name, item in zip(names_list, input_files):
            vcf = pysam.VariantFile(item, 'r')
            rows = 0
            for r in vcf.fetch():
                gk, anchor = _record_group_key(r)
                bi = _blob_index_for(group_order, boundaries, gk, anchor)
                line = str(r)  # full VCF record line incl. trailing newline; columns 0-9 only (one sample)
                if not line.endswith("\n"):
                    line += "\n"
                handles[bi].write(name + "\t" + line)
                rows += 1
            vcf.close()
            input_rows[name] = rows
            if show_progress:
                logging.info(f"Finished splitting {item}")
    finally:
        for h in handles:
            h.close()
    return blob_paths, input_rows


def _open_vcf_text(path):
    """Open a text VCF (plain or gzipped), consume its header, and return the
    file handle positioned at the first data line together with the list of
    sample names declared in the #CHROM header line.
    """
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path, 'r')
    samples = None
    for line in f:
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
            break
    if samples is None:
        f.close()
        return None, None
    return f, samples


def sort_into_single_file(out_f, vcf_header, file_paths_to_combine, sample_list):
    outf = VcfWriter(out_f, target_header=vcf_header)
    file_iterators = []
    file_samples = []
    for item in file_paths_to_combine:
        f, samples = _open_vcf_text(item)
        if f is None:
            continue
        file_iterators.append(f)
        file_samples.append(samples)
    if not file_iterators:
        return 0

    var_q = []
    done_count = 0
    for idx, f in enumerate(file_iterators):
        try:
            line = next(f)
            parts = line.rstrip('\n').split('\t')
            heappush(var_q, (parts[0], int(parts[1]), idx, line))
        except StopIteration:
            file_iterators[idx].close()
            file_iterators[idx] = None
            done_count += 1

    written = 0
    while done_count < len(file_iterators):
        if var_q:
            chrom, pos, idx, line = heappop(var_q)
            # ensure correct sample ordering
            parts = line.rstrip('\n').split('\t')
            main_record = parts[:9]
            samp_records = parts[9:]
            record_samples = {k: v for k, v in zip(file_samples[idx], samp_records)}
            n_fmt = len(parts[8].split(':'))
            rd = []
            c = 0
            for samp in sample_list:
                if samp in record_samples:
                    rd.append(record_samples[samp])
                    c += 1
                else:
                    rd.append('0/0:' + ':'.join(['0'] * (n_fmt - 1)))
            if c < len(record_samples):
                raise RuntimeError('Number of samples in record was greater than out file, please report this')
            main_record += rd
            outf.write('\t'.join(main_record) + '\n')
            written += 1

            if file_iterators[idx] is None:
                continue
            try:
                line = next(file_iterators[idx])
                parts = line.rstrip('\n').split('\t')
                heappush(var_q, (parts[0], int(parts[1]), idx, line))
            except StopIteration:
                file_iterators[idx].close()
                file_iterators[idx] = None
                done_count += 1
        else:
            for idx, f in enumerate(file_iterators):
                if f is None:
                    continue
                try:
                    line = next(f)
                    parts = line.rstrip('\n').split('\t')
                    heappush(var_q, (parts[0], int(parts[1]), idx, line))
                except StopIteration:
                    file_iterators[idx].close()
                    file_iterators[idx] = None
                    done_count += 1
    outf.close()
    assert all(i is None for i in file_iterators)
    return written


def shard_data(args, input_files, Global, show_progress):
    logging.info(f"Merge distance: {args['merge_dist']} bp")
    out_f = args['svs_out'] if args['svs_out'] else '-'
    logging.info("SVs output to {}".format(out_f if out_f != '-' else 'stdout'))

    to_delete = []
    seen_names, names_list = get_names_list(input_files, ignore_csv=False)

    target_rows = args.get("target_rows_per_job", 50000)

    # Pass 1 (parallel): scan every input for its sorted grouping keys, validating sort order.
    if show_progress:
        logging.info("Scanning input files")
    scan_args = [(item, name) for name, item in zip(names_list, input_files)]
    scan_labels = [f"scan sample {name} ({item})" for name, item in zip(names_list, input_files)]
    scan_results = run_process_jobs(scan_sample_keys, scan_args, args['procs'], scan_labels)
    # scan_results come back in job order, which matches names_list/input_files order.
    per_sample_keys = [keys for _, keys, _ in scan_results]

    # Compute blob boundaries: contiguous, ~target_rows-sized spans cut only at safe boundaries.
    group_order, boundaries = compute_blob_boundaries(per_sample_keys, args["merge_dist"], target_rows)
    logging.info(f"Merging into {len(boundaries)} blobs")

    # Capture a reference header (and contig lines) from the first input for output writing.
    ref_header, ref_contigs, _ = _read_vcf_header(input_files[0])
    contig_names = "\n" + ref_contigs.strip() if ref_contigs else ref_contigs

    # Pass 2 (serial): stream all inputs and write each record to its blob file.
    if show_progress:
        logging.info("Writing blobs")
    blob_paths, input_rows = write_blobs(args['wd'], input_files, names_list, group_order, boundaries,
                                         show_progress)
    to_delete += blob_paths

    # Pass 3 (parallel): merge each blob independently.
    if show_progress:
        logging.info("Processing blobs")
    needed_args = {k: args[k] for k in ["wd", "metrics", "merge_within", "merge_dist", "collapse_nearby",
                                        "merge_across", "out_format", "separate", "verbosity", "add_kind",
                                        "progressive"]}
    job_args2 = []
    job_sizes = []
    job_labels2 = []
    merged_outputs = []
    for bi, blob_path in enumerate(blob_paths):
        if os.stat(blob_path).st_size == 0:
            continue
        target_args = needed_args.copy()
        fout = os.path.join(args['wd'], f'blob_{bi}_merged.vcf')
        target_args["svs_out"] = fout
        merged_outputs.append(fout)
        to_delete.append(fout)
        # log_messages=False and show_progress=False: per-job start/finish lines are emitted by
        # run_process_jobs so they stay one-per-job instead of interleaving across worker processes.
        job_args2.append((target_args, [blob_path], seen_names, names_list, False, False,
                          f"blob_{bi}", ref_header, contig_names))
        job_sizes.append(os.path.getsize(blob_path))
        job_labels2.append(f"merge blob_{bi}")

    # Biggest jobs first so a large blob never tails the schedule alone.
    order = sorted(range(len(job_args2)), key=lambda i: job_sizes[i], reverse=True)
    job_args2 = [job_args2[i] for i in order]
    job_labels2 = [job_labels2[i] for i in order]
    if args['procs'] > 1:
        run_process_jobs(process_file_list, job_args2, args['procs'], job_labels2, progress=show_progress)
    else:
        for ji, (ja, label) in enumerate(zip(job_args2, job_labels2)):
            if show_progress:
                logging.info(f"[{ji + 1}/{len(job_args2)}] started {label}")
            process_file_list(*ja)
            if show_progress:
                logging.info(f"[{ji + 1}/{len(job_args2)}] finished {label}")

    # Make a header with all the samples in
    vcf_header = None
    for result_file in merged_outputs:
        vcf_header = pysam.VariantFile(result_file, 'r').header
        tmp = list(filter(None, str(vcf_header).split("\n")))
        new_header = pysam.VariantHeader()
        for line in tmp:
            l = str(line)
            if re.search('^#CHROM', l):
                break
            new_header.add_line(l)
        new_header.add_samples(names_list)
        vcf_header = new_header
        break

    sample_list = list(vcf_header.samples)
    logging.info(f"Samples: {sample_list}")
    written = sort_into_single_file(out_f, vcf_header, merged_outputs, sample_list)
    logging.info("Sample rows before merge {}, rows after {}".format([input_rows[k] for k in names_list], written))
    for item in to_delete:
        if os.path.exists(item):
            os.remove(item)


def get_cosine_similarity(r_format_array, candidates, samp, target_keys):
    res = []
    for index, c in candidates:
        try:
            fmt = c.samples[samp]
        except KeyError:
            # dysgu may have added a _INT to the samp name to make it unique
            dup_key = '_'.join(samp.split('_')[:-1])
            if dup_key in c.samples:
                fmt = c.samples[dup_key]
            else:
                raise KeyError(f'Sample {samp} not in {c.samples}')
        vals2 = np.array([fmt[k] if k in fmt and fmt[k] is not None else 0 for k in target_keys]) + 1e-4
        cs = 1 - abs(cosine(r_format_array, vals2))
        if cs > 0.8:
            res.append((index, c, cs))
    return res


def get_variant_key(r):
    if r.info['SVTYPE'] != 'TRA':
        return r.chrom, r.info['SVTYPE']
    return r.chrom, r.info['CHR2']


def find_similar_candidates(current_cohort_file, variant_table, samp):
    tot = 0
    G = nx.Graph()
    cached = {}
    cohort_index = 0
    for r in current_cohort_file.fetch():
        tot += 1
        key = get_variant_key(r)
        if key not in variant_table:
            cohort_index += 1
            continue

        candidates = [i for i in variant_table[key].search_values(r.pos, r.pos+1)]
        if len(candidates):
            # NMB is skipped for older versions of dysgu
            numeric = {k: v if v is not None else 0 for k, v in r.samples[samp].items() if k != "GT" and k != "NMB"}
            target_keys = numeric.keys()
            r_format_array = np.array(list(numeric.values())) + 1e-4
            candidates = get_cosine_similarity(r_format_array, candidates, samp, target_keys)
            for index, c, cs_similarity in candidates:
                cached[(index, r.chrom)] = c
                G.add_edge(cohort_index, (index, r.chrom), weight=cs_similarity)

        cohort_index += 1

    # Find single best u out-edge and single best v out-edge
    matching_edges_u = {}
    matching_edges_v = {}
    for u in G.nodes():
        if isinstance(u, tuple):
            continue
        best = None
        for v in G.neighbors(u):
            if not isinstance(v, tuple):
                continue
            w = G[u][v]["weight"]
            if not best or w > best[1]:
                best = v, w
        if best:
            if best[0] not in matching_edges_v:
                matching_edges_v[best[0]] = (u, best[1])
            elif best[1] >= matching_edges_v[best[0]][1]:
                del matching_edges_u[matching_edges_v[best[0]][0]]
                matching_edges_v[best[0]] = (u, best[1])
            else:
                continue
            if u not in matching_edges_u:
                matching_edges_u[u] = best
            elif best[1] >= matching_edges_u[u][1]:  # update u edge
                del matching_edges_v[matching_edges_u[u][0]]
                matching_edges_u[u] = best

    matches = {}
    for u, (v, w) in matching_edges_u.items():
        if w != 1:
            matches[u] = cached[v]

    return tot, matches


def recreate_header_with_sample(input_vcf, sample_to_keep):
    # Add only the desired samples to the new header
    new_header = pysam.VariantHeader()
    for record in input_vcf.header.records:
        if record.type != 'SAMPLE':
            new_header.add_record(record)
    for sample in input_vcf.header.samples:
        if sample == sample_to_keep:
            new_header.add_sample(sample)
    return new_header


def copy_samp_from_cohort(samp, cohort_path, new_cohort_path):
    cohort = pysam.VariantFile(cohort_path)
    header = recreate_header_with_sample(cohort, samp)
    temp_c = pysam.VariantFile(new_cohort_path, 'w', header=header)
    for record in cohort.fetch():
        new_record = temp_c.new_record()
        new_record.chrom = record.chrom
        new_record.pos = record.pos
        new_record.id = record.id
        new_record.ref = record.ref
        new_record.alts = record.alts
        for key, value in record.info.items():
            new_record.info[key] = value
        for key, value in record.samples[samp].items():
            new_record.samples[samp][key] = value
        temp_c.write(new_record)
    cohort.close()
    temp_c.close()


def split_cohort_into_target_samples(samples, args):
    job_args = []
    for samp, samp_file in samples.items():
        temp_cohort_path = f"{args['wd']}/temp_cohort_{samp}.vcf"
        job_args.append((samp, args["cohort_update"], temp_cohort_path))
    with multiprocessing.Pool(args['procs']) as pool:
        pool.starmap(copy_samp_from_cohort, job_args)
    return dict((i[0], i[2]) for i in job_args)


def update_target_using_matchings(current_cohort_file, current_out_file, matchings, samp):
    # Second pass update matching SVs
    cohort_index = 0
    for r in current_cohort_file.fetch():
        if cohort_index not in matchings:
            current_out_file.write(r)
            cohort_index += 1
            continue
        vr = matchings[cohort_index]
        cohort_index += 1
        updated = False
        ch_fmt = dict(r.samples[samp].items())
        try:
            items = vr.samples[samp]
        except KeyError:
            # dysgu may have added a _INT to the samp name to make it unique
            dup_key = '_'.join(samp.split('_')[:-1])
            if dup_key in vr.samples:
                items = vr.samples[dup_key]
            else:
                raise KeyError(f'Sample {samp} not in cohort samples {vr.samples}')
        for k, v in items.items():
            if k not in ch_fmt or ch_fmt[k] is None:
                continue
            elif isinstance(v, float):
                if int(v * 100) != int(ch_fmt[k] * 100):
                    ch_fmt[k] = v
                    updated = True
            elif v != r.samples[samp][k]:
                ch_fmt[k] = v
                updated = True
        if updated:
            for k, v in ch_fmt.items():
                r.samples[samp][k] = v

        current_out_file.write(r)


def make_updated_sample_level_vcfs(samp, samp_split_path, samp_file_path, updated_file_path):
    variant_table = {}
    samp_file = pysam.VariantFile(samp_file_path)
    for i, r in enumerate(samp_file.fetch()):
        key = get_variant_key(r)
        if key not in variant_table:
            variant_table[key] = IntervalMap(with_data=True)
        variant_table[key].add(r.pos - 500, r.pos + 500, (i, r))

    for k, v in variant_table.items():
        v.build()

    samp_file.close()

    # First pass, find matching SVs between cohort vcf and input vcf
    current_cohort_file = pysam.VariantFile(samp_split_path)
    tot, matchings = find_similar_candidates(current_cohort_file, variant_table, samp)
    assert samp in samp_split_path

    # Update matched records, save to new file
    current_out_file = pysam.VariantFile(updated_file_path, "w",
                                         header=current_cohort_file.header)
    update_target_using_matchings(current_cohort_file, current_out_file, matchings, samp)
    current_cohort_file.close()
    current_out_file.close()
    os.remove(samp_split_path)


def update_cohort_only(args):

    cohort = pysam.VariantFile(args["cohort_update"])
    cohort_samples = set(cohort.header.samples)
    header = cohort.header
    input_command = ' '.join(sys.argv)
    header.add_line(f'##command="{input_command}"')
    if args["svs_out"] == "-" or args["svs_out"] is None:
        args["svs_out"] = "-"
    outfile = pysam.VariantFile(args["svs_out"] , "w", header=header)

    samples = {}
    samps_counts = defaultdict(int)
    for pth in args["input_files"]:
        v = pysam.VariantFile(pth)
        samps = list(v.header.samples)
        if len(samps) > 1:
            raise ValueError(f"Only one sample supported per input file but {pth} has {list(samps)}")
        samp = samps[0]
        samps_counts[samp] += 1
        if samp in samples:
            logging.warning(f"{samp} already in samples list, {samp} is assumed to be {samp}_{samps_counts[samp]} in cohort file")
            samp = f"{samp}_{samps_counts[samp]}"
        if samp not in cohort_samples:
            raise ValueError(f"Input file has sample name {samp}, but this was not found in the cohort file {args['cohort_update']}")
        samples[samp] = pth

    logging.info("Separating target samples from cohort")
    samp_split_paths = split_cohort_into_target_samples(samples, args)

    logging.info("Processing input samples")
    job_args = []
    for samp, samp_file_path in samples.items():
        job_args.append((samp, samp_split_paths[samp], samp_file_path, f"{args['wd']}/temp_updated_sample_{samp}.vcf"))

    with multiprocessing.Pool(args['procs']) as pool:
        pool.starmap(make_updated_sample_level_vcfs, job_args)

    # Finally write all the cohort shards to output file
    updated_vcfs = {samp: pysam.VariantFile(f"{args['wd']}/temp_updated_sample_{samp}.vcf") for samp in samples.keys()}
    if args["svs_out"] == "-":
        logging.info("SVs output to stdout")
    else:
        logging.info("SVs output to {}".format(args["svs_out"]))

    for cohort_r in cohort:
        updated = False
        for samp, vcf in updated_vcfs.items():
            samp_r = next(vcf)
            cohort_samp = cohort_r.samples[samp]
            for k, v in samp_r.samples[samp].items():
                if v != cohort_samp[k]:
                    updated = True
                    cohort_samp[k] = v
        if updated:
            probs = []
            for s in cohort_samples:
                probs.append(cohort_r.samples[s]["PROB"])
            mean = sum(probs) / len(probs)
            max_prob = max(probs)
            cohort_r.info["MeanPROB"] = mean
            cohort_r.info["MaxPROB"] = max_prob
        outfile.write(cohort_r)

    outfile.close()
    cohort.close()

    for samp in samples.keys():
        tmp_file = f"{args['wd']}/temp_updated_sample_{samp}.vcf"
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

    return


def view_file(args):

    t0 = time.time()

    if args["input_list"]:
        with open(args["input_list"], "r") as il:
            args["input_files"] += tuple([i.strip() for i in il.readlines()])

    if args["separate"] == "True":
        if args["out_format"] == "vcf":
            raise ValueError("--separate only supported for --out-format csv")

        if not all(os.path.splitext(i)[1] == ".csv" for i in args["input_files"]):
            raise ValueError("All input files must have .csv extension")

    if args["cohort_update"]:
        for opt, v in (("out_format", "vcf"), ("merge_within", "False"), ("merge_across", "True"), ("add_kind", "False"), ("separate", "False")):
            if args[opt] != v:
                raise ValueError(f"{opt}={v} fot supported with --cohort-update")
        if args["wd"] is None:
            raise ValueError("Need a working directory --wd")
        update_cohort_only(args)

    else:
        args["metrics"] = False  # only option supported so far
        args["contigs"] = False

        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        manager = multiprocessing.Manager()
        Global = manager.Namespace()
        Global.open_files = soft  # we dont know how many file descriptors are already open by the user, assume all
        Global.soft = soft

        args['progressive'] = args['merge_method'] == "progressive" or (args['merge_method'] == "auto" and len(args['input_files']) > 4)
        prog = '' if args['merge_method'] != "auto" else (': progressive' if args['progressive'] else ': all-vs-all')
        logging.info(f"Merge method {args['merge_method']}{prog}")
        if not args["wd"]:
            if args['procs'] > 1:
                logging.warning(f"A working directory is needed to use multiprocessing with --procs={args['procs']}")
            seen_names, names_list = get_names_list(args["input_files"])
            process_file_list(args, args["input_files"], seen_names, names_list, log_messages=True)
        else:
            shard_data(args, args["input_files"], Global=Global, show_progress=args['progress'])
            if args['clean']:
                os.rmdir(args['wd'])

    logging.info("dysgu merge complete h:m:s, {}".format(str(datetime.timedelta(seconds=int(time.time() - t0))),
                                                    time.time() - t0))
