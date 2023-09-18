import glob
import os
import sys
import pandas as pd
import sortedcontainers
import logging
import time
import datetime
from collections import defaultdict
from dysgu import io_funcs, cluster
from dysgu.map_set_utils import echo, merge_intervals
import multiprocessing
import gzip
import pysam
import numpy as np
from sys import stdout
from heapq import heappop, heappush
import resource
from copy import deepcopy
import re


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


def merge_df(df, n_samples, merge_dist, tree=None, merge_within_sample=False, aggressive=False, log_messages=True):
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
        found = cluster.merge_events(potential, merge_dist, tree, try_rev=False, pick_best=False, add_partners=True,
                                     aggressive_ins_merge=True,
                                     same_sample=False)
        ff = defaultdict(set)
        for f in found:
            if f.partners is None:
                ff[f.event_id] = set([])
            else:
                # Remove partners from same sample
                current = f["table_name"]
                targets = set([])
                passed = True
                for item in f["partners"]:  # Only merge with one row per sample
                    t_name = df.loc[item]["table_name"]
                    if t_name != current and t_name not in targets and len(targets) < n_samples:
                        targets.add(t_name)

                    elif not aggressive:
                        # Merged with self event. Can happen with clusters of SVs with small spacing
                        # e.g. a is merged with b and c, where a is from sample1 and b and c are from sample2
                        # safer not to merge? otherwise variants can be lost
                        passed = False
                if passed:  # enumerate support between components
                    g = f["partners"] + [f["event_id"]]
                    for t1 in g:
                        for t2 in g:
                            if t2 != t1:
                                ff[t1].add(t2)

        df["partners"] = [ff[i] if i in ff else set([]) for i in df.index]
        return df
    else:
        found = cluster.merge_events(potential, merge_dist, tree, try_rev=False, pick_best=True, add_partners=False,
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
            if len(r["partners"]) == 0:
                p2.append("")
                continue
            else:
                key = "|".join([f"{df.loc[idx]['table_name']},{idx}" for idx in r["partners"]])
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
                    ori.append("|".join([f"{df.loc[i]['table_name']},{df.loc[i]['event_id']}" for i in item]))
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


def vcf_to_df(path):
    if os.stat(path).st_size == 0:
        logging.warning(f"File is empty: {path}")
        return pd.DataFrame()

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
    if not header:
        logging.critical(f"File does not have vcf header: {path}")
        quit()

    sample = last.strip().split("\t")[9]
    if path == '-':
        path_h = sys.stdin
    else:
        path_h = path
    df = pd.read_csv(path_h, index_col=None, comment="#", sep="\t", header=None)
    if len(df.columns) > 10:
        msg = f"Can only merge files with one sample in. N samples = {len(df.columns) - 9}"
        raise ValueError(msg)
    parsed = pd.DataFrame()
    parsed["chrA"] = df[0]
    parsed["posA"] = df[1]
    parsed["event_id"] = df[2]
    parsed["ref_seq"] = df[3]
    parsed["variant_seq"] = df[4]
    parsed["filter"] = df[6]
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
               "NMP": ("NMpri", float),
               "NMB": ("NMbase", float),
               "NMS": ("NMsupp", float),
               "MAPQP": ("MAPQpri", float),
               "MAPQS": ("MAPQsupp", float),
               "NP": ("NP", np.int64),
               "OL": ("query_overlap", np.int64),
               "MAS": ("maxASsupp", np.int64),
               "SU": ("su", np.int64),
               "WR": ("spanning", np.int64),
               "PE": ("pe", np.int64),
               "SR": ("supp", np.int64),
               "SC": ("sc", np.int64),
               "BND": ("bnd", np.int64),
               "SQC": ("sqc", float),
               "SCW": ("scw", float),
               "SQR": ("clip_qual_ratio", float),
               "RT": ("type", str),
               "BE": ("block_edge", np.int64),
               "COV": ("raw_reads_10kb", float),
               "MCOV": ("mcov", float),
               "LNK": ("linked", np.int64),
               "CONTIGA": ("contigA", str),
               "CONTIGB": ("contigB", str),
               "ref_seq": ("ref_seq", str),
               "variant_seq": ("variant_seq", str),
               "GC": ("gc", float),
               "NEIGH": ("neigh", np.int64),
               "NEIGH10": ("neigh10kb", np.int64),
               "REP": ("rep", float),
               "REPSC": ("rep_sc", float),
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
               "SBT": ("strand_binom_t", float),
               "PROB": ("prob", float),
               "NG": ("n_gaps", float),
               "NSA": ("n_sa", float),
               "NXA": ("n_xa", float),
               "NMU": ("n_unmapped_mates", np.int64),
               "NDC": ("double_clips", np.int64),
               "RMS": ("remap_score", np.int64),
               "RED": ("remap_ed", np.int64),
               "BCC": ("bad_clip_count", np.int64),
               "STL": ("n_small_tlen", np.int64),
               "RAS": ("ras", np.int64),
               "FAS": ("fas", np.int64),
               "ICN": ("inner_cn", float),
               "OCN": ("outer_cn", float),
               "CMP": ("compress", float),
               "FCC": ("fcc", float),
               "RR": ("ref_rep", float),
               "JIT": ("jitter", float),
               "LEFT_SVINSSEQ": ("left_ins_seq", str),
               "RIGHT_SVINSSEQ": ("right_ins_seq", str),
               }

    df.rename(columns={k: v[0] for k, v in col_map.items()}, inplace=True)
    for value_name, value_type in col_map.values():
        if value_name != "posB_tra" and value_name not in df:
            if value_type == np.int64 or value_type == float:
                df[value_name] = [0] * len(df)
            else:
                df[value_name] = [''] * len(df)

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
                    df[k] = df[k].fillna(0)
                try:
                    df[k] = df[k].astype(dtype)
                except ValueError or OverflowError:
                    raise ValueError("Problem for feature {}, could not interpret as {}".format(k, dtype))
    if "contigA" not in df:
        df["contigA"] = [""] * len(df)
    if "contigB" not in df:
        df["contigB"] = [""] * len(df)
    if 'posB_tra' in df:
        df["posB"] = [i if svt != 'TRA' else j for i, j, svt in zip(df['posB'], df['posB_tra'], df['svtype'])]
        del df['posB_tra']
    df["posA"] = df["posA"].astype(int) - 1  # convert to 0-indexing
    df["posB"] = df["posB"].astype(int) - 1
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


def process_file_list(args, file_list, seen_names, names_list, log_messages):
    dfs = []
    header = None
    contig_names = None
    for bname, item in zip(names_list, file_list):
        name, ext = os.path.splitext(item)

        if ext != ".csv":  # assume vcf
            df, header, _, contig_names = vcf_to_df(item)  # header here, assume all input has same number of fields
        else:
            df = pd.read_csv(item, index_col=None)

        name = list(set(df["sample"]))
        if len(name) > 1:
            msg = f"More than one sample in {item}"
            raise ValueError(msg)

        df["sample"] = [bname] * len(df)
        df["table_name"] = [bname for _ in range(len(df))]
        if len(set(df["sample"])) > 1:
            raise ValueError("More than one sample per input file")

        df = mung_df(df, args)
        if args["merge_within"] == "True":
            l_before = len(df)
            df = merge_df(df, 1, args["merge_dist"], {}, merge_within_sample=True,
                          aggressive=args['collapse_nearby'] == "True", log_messages=log_messages)
            if log_messages:
                logging.info("{} rows before merge-within {}, rows after {}".format(name, l_before, len(df)))
        if "partners" in df:
            del df["partners"]
        dfs.append(df)

    df = pd.concat(dfs)
    if args["merge_across"] == "True":
        if len(seen_names) > 1:
            df = merge_df(df, len(seen_names), args["merge_dist"], {}, aggressive=args['collapse_nearby'] == "True",
                          log_messages=log_messages)

    df = df.sort_values(["chrA", "posA", "chrB", "posB"])

    outfile = open_outfile(args, names_list, log_messages=log_messages)
    if args["out_format"] == "vcf":
        count = io_funcs.to_vcf(df, args, seen_names, outfile, header=header, small_output_f=True,
                                contig_names=contig_names, show_names=log_messages)
        if log_messages:
            logging.info("Sample rows before merge {}, rows after {}".format(list(map(len, dfs)), count))
    else:
        to_csv(df, args, outfile, small_output=False)


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


def shard_job(wd, item_path, name, Global):
    shards = {}
    vcf = pysam.VariantFile(item_path, 'r')
    contigs = len(vcf.header.contigs)
    if contigs == 0:
        contigs = 250
    ub = contigs * contigs
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    ub += soft
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, ub), hard))
    Global.soft = ub

    rows = 0
    for r in vcf.fetch():
        svtype = r.info['SVTYPE']
        svtype = 'INSDUP' if (svtype == "INS" or svtype == "DUP") else svtype
        chrom = r.chrom
        if svtype != "TRA":
            key = f'{svtype}_{chrom}', name
        else:
            chrom2 = r.info["CHR2"]
            key = f'{svtype}_{min(chrom, chrom2)}_{max(chrom, chrom2)}', name
        if key not in shards:
            shards[key] = VcfWriter(os.path.join(wd, name + "~" + key[0] + f".vcf"), vcf.header, name)
        shards[key].write(r)
        rows += 1
    for v in shards.values():
        v.close()
    vcf.close()
    return list(shards.keys()), rows, name


def sort_into_single_file(out_f, vcf_header, file_paths_to_combine, sample_list):
    outf = VcfWriter(out_f, target_header=vcf_header)
    file_iterators = []
    count = 0
    for item in file_paths_to_combine:
        file_iterators.append(pysam.VariantFile(item, 'r'))
        count += 1
    if not file_iterators:
        return 0

    var_q = []
    done_count = 0
    for idx, item in enumerate(file_iterators):
        try:
            v = item.__next__()
            heappush(var_q, (v.chrom, v.pos, idx, v))
        except StopIteration:
            file_iterators[idx].close()
            file_iterators[idx] = None
            done_count += 1

    written = 0
    while done_count < len(file_iterators):
        if var_q:
            chrom, pos, idx, v = heappop(var_q)
            # ensure correct sample ordering
            str_v = str(v).strip().split('\t')
            main_record = str_v[:9]
            samp_records = str_v[9:]
            record_samples = {k: v for k, v in zip(v.samples.keys(), samp_records)}
            rd = []
            c = 0
            for samp in sample_list:
                if samp in record_samples:
                    rd.append(record_samples[samp])
                    c += 1
                else:
                    rd.append('0/0:' + ':'.join(list('0' * (len(v.format.keys()) - 1))))
            if c < len(record_samples):
                raise RuntimeError('Number of samples in record was greater than out file, please report this')
            main_record += rd
            outf.write('\t'.join(main_record) + '\n')
            written += 1

            if file_iterators[idx] is None:
                continue
            try:
                v = file_iterators[idx].__next__()
                heappush(var_q, (v.chrom, v.pos, idx, v))
            except StopIteration:
                file_iterators[idx].close()
                file_iterators[idx] = None
                done_count += 1
        else:
            for idx, item in enumerate(file_iterators):
                if item is None:
                    continue
                try:
                    v = item.__next__()
                    heappush(var_q, (v.chrom, v.pos, idx, v))
                except StopIteration:
                    file_iterators[idx].close()
                    file_iterators[idx] = None
                    done_count += 1
    outf.close()
    assert all(i is None for i in file_iterators)
    return written


def shard_data(args, input_files, Global):
    logging.info(f"Merge distance: {args['merge_dist']} bp")
    out_f = args['svs_out'] if args['svs_out'] else '-'
    logging.info("SVs output to {}".format(out_f if out_f != '-' else 'stdout'))

    to_delete = []
    seen_names, names_list = get_names_list(input_files, ignore_csv=False)

    # Split files into shards
    job_args = []
    for name, item in zip(names_list, input_files):
        job_args.append((args['wd'], item, name, Global))
    pool = multiprocessing.Pool(args['procs'])
    results = pool.starmap(shard_job, job_args)
    input_rows = {}
    job_blocks = defaultdict(list)
    for block_keys, rows, sample_name in results:
        input_rows[sample_name] = rows
        for block_id, _ in block_keys:
            job_blocks[block_id].append(sample_name)
            to_delete.append(os.path.join(args['wd'], block_id + '_merged.vcf'))

    # Set upper bound on open files
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, max(Global.soft, soft) * len(input_files)), hard))

    # Process shards
    job_args2 = []
    needed_args = {k: args[k] for k in ["wd", "metrics", "merge_within", "merge_dist", "collapse_nearby", "merge_across", "out_format", "separate", "verbosity", "add_kind"]}
    merged_outputs = []
    for block_id, names in job_blocks.items():
        job_files = glob.glob(os.path.join(args['wd'], '*' + block_id + '.vcf'))
        to_delete += job_files
        srt_keys = {os.path.basename(i).split("~")[0]: i for i in job_files}
        file_targets = tuple(srt_keys[n] for n in names)
        target_args = deepcopy(needed_args)
        fout = os.path.join(args['wd'], block_id + '_merged.vcf')
        merged_outputs.append(fout)
        target_args["svs_out"] = fout
        job_args2.append((target_args, file_targets, seen_names, names, False))

    if args['procs'] > 1:
        pool.starmap(process_file_list, job_args2)
    else:
        for ja in job_args2:
            process_file_list(*ja)

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

    args["metrics"] = False  # only option supported so far
    args["contigs"] = False

    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    manager = multiprocessing.Manager()
    Global = manager.Namespace()
    Global.open_files = soft  # we dont know how many file descriptors are already open by the user, assume all
    Global.soft = soft

    if not args["wd"]:
        seen_names, names_list = get_names_list(args["input_files"])
        process_file_list(args, args["input_files"], seen_names, names_list, log_messages=True)
    else:
        shard_data(args, args["input_files"], Global=Global)
        if args['clean']:
            os.rmdir(args['wd'])

    logging.info("dysgu merge complete h:m:s, {}".format(str(datetime.timedelta(seconds=int(time.time() - t0))),
                                                    time.time() - t0))
