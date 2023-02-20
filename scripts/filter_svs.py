#/usr/bin/env python3
import glob
import random
import pysam
from collections import defaultdict, Counter
import pandas as pd
# import quicksect
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import click
from dysgu.python_api import vcf_to_df
from dysgu.io_funcs import to_vcf, gen_format_fields, make_main_record
import progressbar
import edlib
from skbio.alignment import StripedSmithWaterman
from cigar import Cigar
# from bisect import bisect

def pair_orientation(aln):
    # Calcutumourss the orientation of a read and its pair
    r1 = aln.is_read1
    rev = aln.is_reverse
    mate_rev = aln.mate_is_reverse
    orientation = None
    if r1:
        if rev:
            if mate_rev:
                orientation = "5to5"
            else:
                orientation = "5to3"
        else:
            if mate_rev:
                orientation = "3to5"
            else:
                orientation = "3to3"
    else:
        if rev:
            if mate_rev:
                orientation = "5to5"
            else:
                orientation = "3to5"
        else:
            if mate_rev:
                orientation = "5to3"
            else:
                orientation = "3to3"
    return orientation


def reciprocal_overlap(a, b):
    # Returns the length of overlapping sequence
    ia = (min(a), max(a))
    ib = (min(b), max(b))
    return max(0, min(ia[1], ib[1]) - max(ia[0], ib[0]))


def path_to_sample(p):
    return os.path.splitext(os.path.basename(p))[0]

def get_clip(r):
    ct = r.cigartuples
    if ct == None:
        return "", None
    elif ct[0][0] == 4:
        return r.query_sequence[:ct[0][1]], True
    elif ct[-1][0] == 4:
        return r.query_sequence[-ct[-1][1]:], False
    else:
        return "", None

def get_clip_sided(r, tag):
    ct = r.cigartuples
    if ct == None:
        return "", False, 0
    if tag == 1 and ct[0][0] == 4:
        return r.query_sequence[:ct[0][1]], False, ct[0][1]
    elif tag == 2 and ct[-1][0] == 4:
        return r.query_sequence[-ct[-1][1]:], False, ct[-1][1]
    elif tag == 0:
        ins = []
        idx = 0
        for op, l in ct:
            if op == 5 or op == 2:  # harc clip or deletion
                continue
            elif op != 1:
                idx += l
            else:  # insertion
                if l >= 15:
                    ins.append((r.query_sequence[idx:idx+l], int(l)))
                    idx += l
        if len(ins) == 1:
            return ins[0][0], False, ins[0][1]
        elif len(ins) == 0:
            return "", False, 0
        elif len(ins) > 1:
            return ins, True, 1
    else:
        return "", False, 0

def get_contig_ins_bases(cont):
    i = 1
    if cont[0].islower():
        while cont[i].islower():
            i += 1
        return cont[:i-1], 1, i
    elif cont[-1].islower():
        while cont[-i].islower():
            i += 1
        return cont[-i:], 2, i
    else:
        return "", 3, 0

def sum_total(cig):
    if cig == None:
        return -1
    else:
        return sum([i[0] for i in list(Cigar(cig).items())])

def sum_match(cig):
    if cig == None:
        return -1
    else:
        return sum([i[0] for i in list(Cigar(cig).items()) if i[1] == "="])

def pos_filter(df, normals, tumours):
    var_types = ['<INV>', '<TRA>', '<DEL>', '<DUP>', '<INS>'] 
    before = {}
    print("BEFORE:")
    before["before"] = {}
    for v in var_types:
        tmp_v = df[df["variant_seq"] == v]
        print(f"Total {v}: {len(tmp_v)}, {len(tmp_v[tmp_v['sample'].isin(normals)])}/{len(tmp_v[tmp_v['sample'].isin(tumours)])} (n/t)")
        before["before"][v] = len(tmp_v[tmp_v["sample"].isin(tumours)])
    tmp_v = df[~df["variant_seq"].isin(var_types)]
    print(f"Total sINS: {len(tmp_v)}, {len(tmp_v[tmp_v['sample'].isin(normals)])}/{len(tmp_v[tmp_v['sample'].isin(tumours)])} (n/t)")
    before["before"]["sINS"] = len(tmp_v[tmp_v["sample"].isin(tumours)])
    # duplicates
    df = df.drop_duplicates(subset=["chrA", "posA", "variant_seq"], keep="first") 
    print("AFTER DUPLUICATES:")
    before["duplicates"] = {}
    for v in var_types:
        tmp_v = df[df["variant_seq"] == v]
        print(f"Total {v}: {len(tmp_v)}, {len(tmp_v[tmp_v['sample'].isin(normals)])}/{len(tmp_v[tmp_v['sample'].isin(tumours)])} (n/t)")
        before["duplicates"][v] = len(tmp_v[tmp_v["sample"].isin(tumours)])
    tmp_v = df[~df["variant_seq"].isin(var_types)]
    print(f"Total sINS: {len(tmp_v)}, {len(tmp_v[tmp_v['sample'].isin(normals)])}/{len(tmp_v[tmp_v['sample'].isin(tumours)])} (n/t)")
    before["duplicates"]["sINS"] = len(tmp_v[tmp_v["sample"].isin(tumours)])
    # close
    df["posA10"] = df["posA"]//25
    df["len100"] = df["svlen"]//25
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "posA10"] = df.loc[df["sample"].isin(normals), "posA10"] - 1
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "posA10"] = df.loc[df["sample"].isin(normals), "posA10"] + 2
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "posA10"] = df.loc[df["sample"].isin(normals), "posA10"] - 1
    df.loc[df["sample"].isin(normals), "len100"] = df.loc[df["sample"].isin(normals), "len100"] - 1
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "len100"] = df.loc[df["sample"].isin(normals), "len100"] + 2
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "posA10"] = df.loc[df["sample"].isin(normals), "posA10"] - 1
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "len100"] = df.loc[df["sample"].isin(normals), "len100"] - 2
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "posA10"] = df.loc[df["sample"].isin(normals), "posA10"] + 2
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")
    df.loc[df["sample"].isin(normals), "len100"] = df.loc[df["sample"].isin(normals), "len100"] + 2
    df = df.drop_duplicates(subset=["chrA", "posA10", "len100", "variant_seq"], keep="first")

    df = df.drop(["posA10", "len100"], axis=1)
    print("AFTER CLOSE:")
    before["close"] = {}
    for v in var_types:
        tmp_v = df[df["variant_seq"] == v]
        print(f"Total {v}: {len(tmp_v)}, {len(tmp_v[tmp_v['sample'].isin(normals)])}/{len(tmp_v[tmp_v['sample'].isin(tumours)])} (n/t)")
        before["close"][v] = len(tmp_v[tmp_v["sample"].isin(tumours)])
    tmp_v = df[~df["variant_seq"].isin(var_types)]
    print(f"Total sINS: {len(tmp_v)}, {len(tmp_v[tmp_v['sample'].isin(normals)])}/{len(tmp_v[tmp_v['sample'].isin(tumours)])} (n/t)")
    before["close"]["sINS"] = len(tmp_v[tmp_v["sample"].isin(tumours)])
    print("Removing normals...")
    df = df[df["sample"].isin(tumours)]
    if df.empty:
        print("nope")
        exit()
    return df, before

@click.command()
@click.argument("norm")
@click.argument("tumr")
@click.option("-c", "--csv", default=None)
@click.option("-o", "--outdir")
@click.option("-p", "--pool-glob", default="*.*am")
@click.option("-n", "--num-normals", default=20)
@click.option("-d", "--pairs-df")
@click.option("--normal-col", default="normal")
@click.option("--tumour-col", default="tumour")
@click.option("--normal-col-conv", default=None)
@click.option("--tumour-col-conv", default=None)
@click.option("-i", "--inc-ins", default=True, is_flag=True)
@click.option("--out-name", default="out_all.csv", is_flag=True)
def filter_parental_bams(norm, tumr, csv, outdir, pool_glob, num_normals, pairs_df, normal_col, tumour_col, normal_col_conv, tumour_col_conv, inc_ins, out_name):
    print(norm)
    print(tumr)
    print(outdir)
    print(pool_glob)
    print(pairs_df)
    if os.path.isfile(norm):
        normals = [norm]
        normal_samples = [path_to_sample(norm)]
    if os.path.isfile(tumr):
        tumours = [tumr]
        tumour_samples = [path_to_sample(tumr)]
    if os.path.isdir(norm):
        normals = os.listdir(norm)
        normals = [os.path.join(norm, i) for i in normals]
        normal_samples = [path_to_sample(i) for i in normals]
    if os.path.isdir(tumr):
        tumours = os.listdir(tumr)
        tumours = [os.path.join(tumr, i) for i in tumours]
        tumour_samples = [path_to_sample(i) for i in tumours]
    if pairs_df != None:
        pairs_conv = pd.read_csv(pairs_df)
        pairs_conv = pairs_conv.reset_index()
        print(pairs_conv)
        if normal_col_conv != None:
            normal_samples = [pairs_conv[pairs_conv[normal_col] == i][normal_col_conv].iloc[0] for i in normal_samples]
            print("normal samples: ", normal_samples)
            #normal_col = normal_col_conv
        if tumour_col_conv != None:
            tumour_conv = {}
            for i in tumour_samples:
                tumour_conv[pairs_conv[pairs_conv[tumour_col] == i][tumour_col_conv].iloc[0]] = i
            print(tumour_conv)
            tumour_samples = [pairs_conv[pairs_conv[tumour_col] == i][tumour_col_conv].iloc[0] for i in tumour_samples]
            print("tumour samples: ", tumour_samples)
            #tumour_col = tumour_col_conv
        else:
            tumour_conv = {}
            for i in tumour_samples:
                tumour_conv[i] = i
            print(tumour_conv)
            #tumour_samples = [pairs_conv[pairs_conv[tumour_col] == i][tumour_col_conv].iloc[0] for i in tumour_samples]
            print("tumour samples: ", tumour_samples)
            #tumour_col = tumour_col_conv


    if csv == None:
        print("Reading and combining vcfs...")
        for z, i in enumerate(normals):
            print(f"{round(z/(len(normals)+len(tumours))*100, 2)}%", end="\r")
            if i == normals[0]:
                df = vcf_to_df(i)[0] 
            else:
                n = vcf_to_df(i)
                df = pd.concat([df, n[0]])
        normal_samples = list(set(df["sample"]))
        print("normal samples: ", normal_samples)
        for z, i in enumerate(tumours):
            print(f"{round((z+len(normals))/(len(normals)+len(tumours))*100, 2)}%", end="\r")
            t = vcf_to_df(i)
            df = pd.concat([df, t[0]])
        tumour_samples = list(set(df["sample"]))
        tumour_samples = [i for i in tumour_samples if i not in normal_samples]
        print("tumour samples :", tumour_samples)
        prev_tag = False
        #df.to_csv("shortcut.csv")
    else:
        if csv.endswith(".csv"):
            df = pd.read_csv(csv, sep=",")
        elif csv.endswith(".tsv"):
            df = pd.read_csv(csv, sep="\t")
        else:
            df = pd.read_csv(csv, delim_whitespace=True)
        if "unique" in list(df.columns):
            df = df.rename(columns={"unique": "prev_unique"})
            prev_tag = True
        else:
            prev_tag = False
        df[["contigA", "contigB"]] = df[["contigA", "contigB"]].fillna("")
    print("")
    df = df[df["chrA"].str.contains("hs37d5")==False]
    df = df[df["chrB"].str.contains("hs37d5")==False]
    print("all samples")
    print(df)
    all_samples = set(df["sample"].tolist())
    print(all_samples)
    total = len(df)
    df, before = pos_filter(df, normal_samples, tumour_samples)
    df[["left_ins_seq", "right_ins_seq"]] = df[["left_ins_seq", "right_ins_seq"]].fillna("")
    print("duplicates and normals removed")
    print(df)
    sample_pools = {}
    print(f"Samples: {len(all_samples)}")
    print(all_samples)
    if num_normals == "all":
        num_normals = len(all_samples)
    if pairs_df == None:
        print(f"Using {num_normals} normals")
        for k in all_samples:
            if k in normal_samples:
                sample_pools[k] = [] #[i for i in normals if i != k]
            elif k in tumour_samples:
                sample_pools[k] = random.sample([i for i in normal_samples if i != k], num_normals)
    elif os.path.isfile(pairs_df) and num_normals > 0:
        pairs = pd.read_csv(pairs_df)
        pairs = pairs[[normal_col, tumour_col]]
        for k in all_samples:
            try:
                norm_pair = pairs[pairs[tumour_col]==k][normal_col].values[0]
                tmp_normal_samples = normal_samples.copy()
                tmp_normal_samples.remove(norm_pair)
                if k in normal_samples:
                    sample_pools[k] = [] #[i for i in normals if i != k]
                elif k in tumour_samples:
                    norm_file = os.path.join(os.path.dirname(pool_glob), norm_pair+".bam")
                    if os.path.exists(norm_file):
                        sample_pools[k] = [norm_pair] + random.sample([i for i in tmp_normal_samples if i != k], num_normals-1)
                    else:
                        sample_pools[k] = random.sample([i for i in normal_samples if i != k], num_normals)
            except:
                if k in normal_samples:
                    sample_pools[k] = [] #[i for i in normals if i != k]
                elif k in tumour_samples:
                    sample_pools[k] = random.sample([i for i in normal_samples if i != k], num_normals)

    elif os.path.isfile(pairs_df) and num_normals == 0:
        print("Using pairs")
        pairs = pd.read_csv(pairs_df)
        pairs = pairs[[normal_col, tumour_col]]
        pairs = list(pairs.itertuples(index=False, name=None))
        for k in all_samples:
            for p in pairs:
                if k == p[0]:
                    sample_pools[k] = []
                if k == p[1]:
                    sample_pools[k] = [p[0]]
        print(sample_pools)
    else:
        print(f"{pairs_df} is not a file")
        exit

    if len(normal_samples) == 1 and len(tumour_samples) == 1 and pairs_df != None:
        normal_samples = pairs[normal_col].tolist()
        tumour_samples = pairs[tumour_col].tolist()
        normal_samples = [i.split(".")[0] for i in normal_samples]
        tumour_samples = [i.split(".")[0] for i in tumour_samples]

    for k, v in sample_pools.items():
        print(k, v)
    
        if "*" in pool_glob:
            crams = {path_to_sample(k):
                         k for k in glob.glob(pool_glob)}
        else:
            crams = {path_to_sample(k):
                        k for k in pool_glob.split(",")}
        print(crams)

    cram_pool = {}
    for k, v in crams.items():
        try:
            cram_pool[k] = pysam.AlignmentFile(v, "rb")
        except:
            print(f"{k} {v} failed..")
    unique = []
    c = 0
    bar = progressbar.ProgressBar(maxval=len(df)).start()
    # df.to_csv("pos_filtered_shortcut.csv")
    for count, (idx, r) in enumerate(df.iterrows()):
        if r["sample"] in normal_samples:
            unique.append(None)
            c += 1
            continue
        if prev_tag == True:
            if r["prev_unique"] == True:
                unique.append(False)
                c += 1
                continue
        try:
            pool = sample_pools[r["sample"]]
        except:
            unique.append(1)
            c += 1
            continue 
        pad = (1250, 1250)
        uni = True
        
        chr1, chr2 = r["chrA"], r["chrB"]
        pos1, pos2 = r["posA"], r["posB"]
        var = r["variant_seq"]
        samp = r["sample"]
        #['<INV>', '<TRA>', '<INS>', '<DEL>', '<DUP>']
        not_insertions = {'<INV>', '<TRA>', '<DEL>', '<DUP>'}
        if var in not_insertions:
            if chr1 == chr2 and abs(pos2 - pos1) >= 2500:  # For larger intrachromsomal SVs
                # Fetch the reads from other sample .bam files, just get from one side of breaksite
                for b in pool:
                    try:
                        cram_f = cram_pool[path_to_sample(b)]  # pysam.AlignmentFile(crams[b], "rc")
                    except:
                        continue
                    # Collect alignments from same site in other samples: (chrom1, position1)
                    for aln in cram_f.fetch(chr1, 0 if pos1 - pad[0] < 0 else pos1 - pad[0], pos1 + pad[1]):
                        if all([cram_f.getrname(aln.rnext) == chr2,  # Same chrom
                                not aln.flag & 2,  # Is discordant
                                pos2 - pad[0] < aln.pnext < pos2 + pad[1],  # Pnext is the same
                                reciprocal_overlap((pos1, pos2), (aln.pos, aln.pnext)) > 0,   # Some overlap
                                reciprocal_overlap((pos1, pos2), (aln.pos, aln.pnext)) / float(abs(pos1 - pos2)) > 0.8]):
                                # Overlap is reciprocal
                            uni = False
                            break
                        else:
                            uni = True  # All reads look unique

                    if not uni:
                        break

            
            elif chr1 == chr2 and abs(pos2 - pos1) < 2500:  # For smaller intrachromosomal SVs
                orientation_most_common = r["join_type"] 
                gap = r["svlen"] 
                if gap >= 300:
                    for b in pool:
                        try:
                            cram_f = cram_pool[path_to_sample(b)]
                        except:
                            continue
                        #cram_f = pysam.AlignmentFile(bams[b], "rc")
                        for aln in cram_f.fetch(chr1, 0 if pos1 - 1000 < 0 else pos1 - 1000,
                                               pos1 + 1000):  # Dont collect all reads at once
                            if aln.flag & 2 or cram_f.getrname(aln.rnext) != chr1:
                                continue  # Skip concordant pairs
                            # Could use reciprocal overlap / relative overlap,
                            # but using the pair orientation and the gap size is basically the same
                            if all([abs(aln.pos - aln.pnext) / gap > 0.90,
                                    orientation_most_common == pair_orientation(aln),  # Orientation must be same
                                    abs(aln.pos - aln.pnext) < 3000]):  # Same if insert is 90% of the max size
                                uni = False
                                break

                            else:
                                uni = True

                        if not uni:
                            break
                else:
                    uni = True
                    for b in pool:
                        try:
                            cram_f = cram_pool[path_to_sample(b)]
                        except:
                            continue
                        cigartuples = [aln.cigartuples for aln in cram_f.fetch(chr1, 0 if pos1 - 1000 < 0 else pos1 - 1000,
                                                   pos1 + 1000) if aln.cigartuples != None] 
                        for ctt in cigartuples:
                            seen = any( (0 < ct[0] <= 2 and ct[1] >= 15) or (ct[0] == 4 and ct[1] >= 15) for ct in ctt)
                            if seen == True:
                                uni = False
                                #print(var, ctt)
                                break
                        if not uni:
                            break

            elif chr1 != chr2:  # Deal with translocations
                pad = (2000, 2000)
                if r["su"] >= 3:
                    for b in pool:
                        try:
                            cram_f = cram_pool[path_to_sample(b)]
                        except:
                            continue
                        # Collect alignments from same site in other samples: (chrom1, position1)
                        for aln in cram_f.fetch(chr1, 0 if pos1 - pad[0] < 0 else pos1 - pad[0], pos1 + pad[0]):
                            if cram_f.getrname(aln.rnext) == chr2 and pos2 - pad[0] < aln.pnext < pos2 + pad[1]:
                                uni = False                            
                                break
                            else:
                                uni = True

                        if uni == False:
                            break

        else: 
            if var == "<INS>":# or var.startswith("<") == False:
                conts = []
                for contig in (r["contigA"], r["contigB"]): # r["left_ins_seq"], r["right_ins_seq"]):
                    if len(contig) == 0:
                        continue
                    cont, tag, cont_len = get_contig_ins_bases(contig)
                    if cont_len > 15:
                        conts.append((cont, tag, cont_len))
                if len(conts) == 0:
                    if len(r["left_ins_seq"]) > 0:
                        cont = r["left_ins_seq"]
                        tag = 1
                        cont_len = len(r["left_ins_seq"])
                    if len(r["right_ins_seq"]) > 0:
                        cont = r["right_ins_seq"]
                        tag = 2
                        cont_len = len(r["right_ins_seq"])


                if len(conts) == 0:
                    assert c == count
                    c += 1
                    unique.append(False)
                    bar.update(c)
                    print("lost")
                    continue

                svlen = r["svlen"] 
                lower_bound = pos1-svlen
                upper_bound = pos1+svlen
                for b in pool:
                    try:
                        cram_f = cram_pool[path_to_sample(b)]  # pysam.AlignmentFile(crams[b], "rc")
                    except:
                        continue
                    not_found = True

                    reads = []
                    for aln in cram_f.fetch(chr1, 0 if lower_bound < 0 else lower_bound, upper_bound):
                        reads.append(aln)

                    for aln in reads:
                        for cvs in conts:
                            query = StripedSmithWaterman(cvs[0])
                            clip, multiple_clips, clip_len = get_clip_sided(aln, cvs[1])
                            if multiple_clips == False and clip_len > 0:
                                if clip_len > 10:
                                    al = query(clip)
                                    if al.optimal_alignment_score >= 40:
                                        not_found = False
                                        break
                            elif multiple_clips == True:
                                for clips in clip:
                                    if clips[1] > 10:
                                        al = query(clips[0])
                                        if al.optimal_alignment_score >= 40:
                                            not_found = False
                                            break

                    uni = not_found
                    if not_found == False:
                        break

            else:
                ins_seq = var
                if ins_seq == "":
                    assert c == count
                    c += 1
                    unique.append(False)
                    bar.update(c)
                    continue
                
                not_found = True
                query = StripedSmithWaterman(str(ins_seq))
                for b in pool:
                    try:
                        cram_f = cram_pool[path_to_sample(b)]  # pysam.AlignmentFile(crams[b], "rc")
                    except:
                        continue
                    for aln in cram_f.fetch(chr1, 0 if pos1 - 100 < 0 else pos1 - 100, pos1 + 100):
                        clip, multiple_clips, clip_len = get_clip_sided(aln, 0)
                        if multiple_clips == False and clip_len > 10:
                            alignment = query(str(clip))
                            if alignment.optimal_alignment_score >= 40:
                                not_found = False
                                break
                        if multiple_clips == True: # ~25% faster than yielding clips
                            for clips in clip:
                                if clips[1] > 10: 
                                    alignment = query(str(clips[0]))
                                    #print(alignment)
                                    if alignment.optimal_alignment_score >= 10:
                                        not_found = False
                                        break

                    uni = not_found
                    if not_found == False:
                        break

        assert c == count  # Check no rows were dropped
        c += 1
        unique.append(uni)
        if uni:
            print(r["sample"], r["chrA"], r["posA"], r["chrB"], r["posB"], r["svtype"])

        bar.update(c)
    
    print(df)
    pickle.dump(unique, open("unique.pkl", "wb"))
    df["unique"] = unique
    print(df)
    ba_dict = {}
    var_types = ['<INV>', '<TRA>', '<DEL>', '<DUP>', '<INS>']
    for i in tumour_samples:
        print(i)
        tmp = df[df["sample"] == i]
        ba_dict[i] = (len(tmp), len(tmp[tmp["unique"] == True]))
        print("before: ", len(tmp), " after: ", len(tmp[tmp["unique"] == True]))
        tmp_df = tmp[tmp["unique"] == True] 
        for v in var_types:
            tmp_v = tmp_df[tmp_df["variant_seq"] == v]
            print(f"Total {v}: {len(tmp_v)}")
        tmp_v = tmp_df[~tmp_df["variant_seq"].isin(var_types)]
        print(f"Total sINS: {len(tmp_v)}")

        tmp = tmp[tmp["unique"] == True]
        ins_tmp = tmp[tmp["variant_seq"] == "<INS>"]
        with open(os.path.join(outdir, f"{tumour_conv[i]}.vcf"), "w") as f:
            to_vcf(tmp, {"add_kind": True, "verbosity": 1, "metrics": True}, [tumour_conv[i]], f)

    print("After read based filter:")
    before["pos"] = {}
    tmp_df = df[df["unique"] == True] 
    for v in var_types:
        tmp_v = tmp_df[tmp_df["variant_seq"] == v]
        print(f"Total {v}: {len(tmp_v)} ({before['close'][v]-len(tmp_v)})")
        before["pos"][v] = len(tmp_v)
    tmp_v = tmp_df[~tmp_df["variant_seq"].isin(var_types)]
    print(f"Total sINS: {len(tmp_v)} ({before['close']['sINS']-len(tmp_v)})")
    before["pos"]["sINS"] = len(tmp_v)

    print(ba_dict)
    print([ba_dict[i][0]-ba_dict[i][1] for i in ba_dict])
    print(before)
    stages = ["before", "duplicates", "close", "pos"]
    for v in var_types + ["sINS"]:
        print(f"{v}:")
        for p, i in enumerate(stages):
            if p == 0:
                print(f"{i}: {before[i][v]}")
            else:
                print(f"{i}: {before[i][v]} ({int(before[stages[p-1]][v]) - int(before[i][v])})")
    df.to_csv(os.path.join(outdir, out_name), sep="\t")


if __name__ == "__main__":
    filter_parental_bams()
