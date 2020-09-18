#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True


import numpy as np
cimport numpy as np

from dysgu.map_set_utils cimport unordered_map
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement

from libc.math cimport fabs as c_fabs
from libcpp.vector cimport vector as cpp_vector

from dysgu.map_set_utils import echo

import math
import array
import pickle
import os
import click
import gzip
import pandas as pd
pd.options.mode.chained_assignment = None


class BadClipCounter:

    def __init__(self, n_references):
        self.clip_pos_arr = [array.array("L", []) for i in range(n_references)]  # 32 bit unsigned int

    def add(self, chrom, position):
        self.clip_pos_arr[chrom].append(position)

    def sort_arrays(self):
        self.clip_pos_arr = [np.array(i, dtype=np.dtype("i"))for i in self.clip_pos_arr]
        [i.sort() for i in self.clip_pos_arr]

    def count_near(self, int chrom, int start, int end):

        cdef int [:] a = self.clip_pos_arr[chrom]
        if len(a) == 0:
            return 0

        idx = np.searchsorted(a, start)

        # search forward and backwards
        cdef int i = idx
        cdef int len_a = len(a)
        cdef int count = 0
        cdef int p

        while True:
            if i < 0 or i == len_a:
                break
            p = a[i]
            if p <= end:
                count += 1
                i += 1
            else:
                break

        return count

def get_badclip_metric(events, bad_clip_counter, bam):

    bad_clip_counter.sort_arrays()
    new_events = []
    for e in events:

        if e["chrA"] == e["chrB"] and abs(e["posB"] - e["posA"]) < 500:
            start = min(e["posA"], e["posB"]) - 500
            end = max(e["posA"], e["posB"]) + 500
            count = bad_clip_counter.count_near(bam.gettid(e["chrA"]), start, end)
        else:
            c1 = bad_clip_counter.count_near(bam.gettid(e["chrA"]), e["posA"] - 500, e["posA"] + 500)
            c2 = bad_clip_counter.count_near(bam.gettid(e["chrB"]), e["posB"] - 500, e["posB"] + 500)
            count = c1 + c2
        e["bad_clip_count"] = count
        new_events.append(e)

    return new_events


cdef float soft_clip_qual_corr(reads):
    """Function to compute the correlation of quality values between the soft-clipped portions of reads. values of
    1 imply no correlation whereas closer to 0 implies that quality values are correlated"""
    cdef const unsigned char[:] quals
    cdef int idx, x, i, j
    cdef unordered_map[int, cpp_vector[int]] qq

    for r in reads:

        quals = r.query_qualities
        if r.cigartuples is None or quals is None:
            continue
        if r.cigartuples[0][0] == 4:
            idx = 0
            for x in range(r.pos - r.cigartuples[0][1], r.pos):
                qq[x].push_back(quals[idx])
                idx += 1
        if r.cigartuples[-1][0] == 4:
            end_pos = r.query_alignment_end + r.cigartuples[-1][1]
            for idx in range(-1, - r.cigartuples[-1][1], -1):
                qq[end_pos].push_back(quals[idx])
                end_pos -= 1

    if qq.empty():
        return -1

    cdef float all_z = 0
    cdef bint seen = False
    cdef cpp_vector[int] all_v
    cdef float sum_all_v = 0
    cdef cpp_vector[int] second
    cdef int first
    cdef float mean_second, sum_second, sum_diff
    cdef unordered_map[int, cpp_vector[int]].iterator qq_iter

    with nogil:
        qq_iter = qq.begin()
        while qq_iter != qq.end():

            second = dereference(qq_iter).second

            if second.size() > 1:
                sum_second = 0
                for i in range(second.size()):
                    sum_second += second[i]
                mean_second = sum_second / float(second.size())
                for i in range(second.size()):
                    all_z += c_fabs(second[i] - mean_second)

                sum_all_v += sum_second
                all_v.insert(all_v.end(), second.begin(), second.end())  # concat

                seen = True
            preincrement(qq_iter)

    if not seen:
        return -1

    cdef float mean_all_v = sum_all_v / float(all_v.size())

    cdef float z = 0
    for i in range(all_v.size()):
        z += c_fabs(all_v[i] - mean_all_v)

    if z == 0:
        return -1

    return all_z / z


# stolen from svtyper with minor modifications
# https://github.com/hall-lab/svtyper/blob/master/svtyper/singlesample.py
# efficient combinatorial function to handle extremely large numbers
def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in range(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# return the genotype and log10 p-value
def bayes_gt(ref, alt, is_dup):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [1e-2, 0.2, 1/3.0]
    else:
        p_alt = [1e-3, 0.5, 0.9]

    total = ref + alt
    log_combo = log_choose(total, alt)

    lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return lp_homref, lp_het, lp_homalt


def get_gt_metric(events, infile, add_gt=False):
    if add_gt:
        click.echo("Adding GT", err=True)
    new_events = []
    pad = 50
    for e in events:

        if infile is None or not add_gt:
            e["GT"] = "./."
            e["GQ"] = -1
            continue

        regions = [[e["chrA"], e["posA"] - pad, e["posA"] + pad]]

        if e["svtype"] != "INS":
            sep = abs(e["posB"] - e["posA"])
            if e["chrB"] == e["chrA"] and sep < 1000:
                if e["posB"] >= e["posA"]:
                    regions[0][-1] = e["posB"] + pad
                else:
                    regions[0][0] = e["posB"] - pad
            else:
                regions.append([e["chrB"], e["posB"] - pad, e["posB"] + pad])

        d = set([])
        probable_sv_reads = set([])

        svlen = e["svlen"] + 1e-6
        wr = e["spanning"] > 0
        deletion = e["svtype"] == "DEL"
        posA = []

        fail = False
        for chrom, start, end in regions:

            target_pos = start + pad
            try:
                fetch = infile.fetch(chrom, start, end)
            except:
                fail = True
                break

            for r in fetch:
                ct = r.cigartuples

                if ct and r.flag & 2 and r.mapq > 10:
                    if not r.pos <= target_pos <= r.reference_end:
                        continue
                    if wr or e["svlen"] < 50:
                        current_pos = r.pos
                        first = True
                        for opp, length in ct:
                            if opp == 4 or opp == 5:
                                if (first and abs(target_pos - current_pos) < 100) or \
                                    (not first and abs(current_pos + length - target_pos) < 100):
                                    probable_sv_reads.add(r.qname)
                                    break
                            elif opp == 2:
                                if deletion and abs(target_pos - current_pos) < 100 and min(svlen, length) / max(svlen, length) > 0.8:
                                    probable_sv_reads.add(r.qname)
                                    break
                                current_pos += length
                            elif opp == 1:
                                if not deletion and abs(target_pos - current_pos) < 100 and min(svlen, length) / max(svlen, length) > 0.8:
                                    probable_sv_reads.add(r.qname)
                                    break
                                current_pos += 1

                            elif opp == 0 or opp == 7 or opp == 8 or opp == 3:  # All match, match (=), mis-match (X), N's
                                current_pos += length
                            first = False
                        else:
                            d.add(r.qname)
                    else:
                        first, last = ct[0][0], ct[-1][0]
                        if first == 4 or first == 5 or last == 4 or last == 5:
                            probable_sv_reads.add(r.qname)
                        else:
                            d.add(r.qname)

        if fail or len(d) == 0 and len(probable_sv_reads) == 0:
            e["GT"] = "./."
            e["GQ"] = -1
            continue

        ref = len(d.difference(probable_sv_reads))

        support_reads = e["su"] - e["spanning"]  # spanning are counted twice

        gt_lplist = bayes_gt(ref, support_reads, e["svtype"] == "DUP")
        best, second_best = sorted([ (i, e) for i, e in enumerate(gt_lplist) ], key=lambda x: x[1], reverse=True)[0:2]
        gt_idx = best[0]

        gt_sum = 0
        result = {}
        for gt in gt_lplist:
            try:
                gt_sum += 10**gt
            except OverflowError:
                gt_sum += 0
        if gt_sum > 0:
            gt_sum_log = math.log(gt_sum, 10)
            sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
            phred_gq = min(-10 * (second_best[1] - best[1]), 200)
            result['GQ'] = int(phred_gq)
            result['SQ'] = sample_qual
            if gt_idx == 1:
                result['GT'] = '0/1'
            elif gt_idx == 2:
                result['GT'] = '1/1'
            elif gt_idx == 0:
                result['GT'] = '0/0'
        else:
            result['GQ'] = '.'
            result['SQ'] = '.'
            result['GT'] = './.'
        e.update(result)

    return events


def apply_model(df, mode, contigs, paired, thresholds):

    pth = os.path.dirname(os.path.abspath(__file__))
    pth = f"{pth}/dysgu_model.1.pkl.gz"
    models = pickle.load(gzip.open(pth, "rb"))

    if paired == "False":
        col_key = f"pacbio_cols"
        model_key = f"pacbio_classifier"
        if contigs == "False":
            col_key += "_no_contigs"
            model_key += "_no_contigs"
    else:
        col_key = f"pe_cols"
        model_key = f"pe_classifier"

    cols = models[col_key]
    clf = models[model_key]

    c = dict(zip(
        ['SQC', 'CIPOS95',  'GC', 'REP', 'REPSC',  'SU', 'WR',       'SR',   'SC', 'NEXP',        'RPOLY',          'STRIDE', 'SVTYPE', 'SVLEN', 'NMP',   'NMB',    'MAPQP',   'MAPQS',    'NP', 'MAS',       'BE',         'COV',            'MCOV', 'NEIGH', 'NEIGH10',   'RB',        'PS',   'MS',    'NG',     'NSA',  'NXA',  'NMU',              'NDC',          'RMS',         'RED',      'BCC',            'STL'],
        ['sqc', 'cipos95A', 'gc', 'rep', 'rep_sc', 'su', 'spanning', 'supp', 'sc', 'n_expansion', 'ref_poly_bases', 'stride', 'svtype', 'svlen', 'NMpri', 'NMbase', 'MAPQpri', 'MAPQsupp', 'NP', 'maxASsupp', 'block_edge', 'raw_reads_10kb', 'mcov', 'neigh', 'neigh10kb', 'ref_bases', 'plus', 'minus', 'n_gaps', 'n_sa', 'n_xa', 'n_unmapped_mates', 'double_clips', 'remap_score', 'remap_ed', 'bad_clip_count', 'n_small_tlen']
                 ))
    echo(c)
    X = df[[c[i] for i in cols]]
    X.columns = cols
    keys = {"DEL": 1, "INS": 2, "DUP": 3, "INV": 4, "TRA": 2, "INV:DUP": 2}

    X["SVTYPE"] = [keys[i] for i in X["SVTYPE"]]
    X["SVLEN"] = [i if i == i else -1 for i in X["SVLEN"]]

    for c in models["cats"]:  # categorical data
        if c in X:
            X[c] = [i if i == i else 0 for i in X[c]]
            X[c] = X[c].astype("category")

    pred = np.round(models[model_key].predict_proba(X)[:, 1], 3)
    df = df.assign(prob=pred)

    df = df.assign(filter=["PASS" if ((svt in thresholds and i >= thresholds[svt]) or (svt not in thresholds and i >= 0.5)) else "lowProb" for svt, i in zip(df["svtype"], pred)])

    return df
