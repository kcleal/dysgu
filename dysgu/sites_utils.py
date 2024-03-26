from collections import namedtuple, deque, Counter
import pysam
import pandas as pd
import logging
from dysgu.map_set_utils import echo


def parse_variant_seqs_dysgu(r, svt, chrom, start, chrom2, stop, paired_end, ref_genome):

    a_left = None
    a_right = None

    if not paired_end:  # long-read analysis doesnt use sequence clustering
        return a_left, a_right

    if stop == start:
        stop += 1
    if svt == "INS":
        alt = r.alts[0]
        if alt[0] != "<":
            a_left = alt[0]  # for clustering left-hand-side soft clipped reads
            a_right = a_left

        elif "LEFT_SVINSSEQ" in r.info or "RIGHT_SVINSSEQ" in r.info:
            if "LEFT_SVINSSEQ" in r.info and len(r.info["LEFT_SVINSSEQ"]) > 0:
                a_left = r.info["LEFT_SVINSSEQ"]
            if "RIGHT_SVINSSEQ" in r.info and len(r.info["RIGHT_SVINSSEQ"]) > 0:
                a_left = r.info["RIGHT_SVINSSEQ"]

    elif svt == "DEL":
        if len(r.ref) > 1:
            a_left = r.ref
            a_right = a_left
        else:
            if stop - start < 100:

                a_left = ref_genome.fetch(chrom, start, stop)
                a_right = a_left
            else:
                a_left = ref_genome.fetch(chrom, start, start + 100)
                a_right = ref_genome.fetch(chrom, stop - 100, stop)

    return a_left, a_right


Site = namedtuple("Site", ["chrom", "start", "chrom2", "end", "svtype", "index", "id", "svlen", "prob"])


def vcf_reader(pth, infile, parse_probs, sample_name, ignore_sample, default_prob=0.6, pass_only=True):
    if pth is None:
        return None

    logging.info(f"Reading --sites {pth}")
    vcf = pysam.VariantFile(pth)

    samples = list(vcf.header.samples)
    if len(samples) == 1 or sample_name not in samples:
        ignore_sample = False

    # variants must be fed into graph in roughly sorted order
    recs = {}
    not_parsed = []
    for idx, r in enumerate(vcf):
        if pass_only and "PASS" not in r.filter:
            continue
        if ignore_sample and sample_name in r.samples and "GT" in r.samples[sample_name]:
            this_gt = r.samples[sample_name]['GT']
            if not (this_gt == (0, 0) or this_gt == "0/0" or this_gt == "0" or this_gt == "." or this_gt == "./." or this_gt == "0|0"):
                continue
        if "CHROM2" in r.info:
            chrom2_info = r.info["CHROM2"]
        else:
            chrom2_info = r.chrom
        svt = r.info["SVTYPE"]
        if r.chrom != chrom2_info and svt in {"INS", "DEL", "DUP", "INV"}:
            raise ValueError(f"CHROM2 must equal chrom for SVTYPE {svt}, (chrom {r.chrom}, pos {r.start})")
        if r.chrom != chrom2_info and svt == "BND":
            svt = "TRA"

        if "DUP" in svt and svt != "DUP":
            svt = "DUP"

        if svt not in {"INS", "DEL", "DUP", "INV", "TRA", "BND"}:
            not_parsed.append(svt)
            continue

        chrom = infile.gettid(r.chrom)

        # try and switch between "chr" representation
        if chrom == -1:
            if "chr" in r.chrom:
                chrom = infile.gettid(r.chrom[3:])
            else:
                chrom = infile.gettid("chr" + r.chrom)

        chrom2 = infile.gettid(chrom2_info)
        if chrom2 == -1:
            if "chr" in chrom2_info:
                chrom2 = infile.gettid(chrom2_info[3:])
            else:
                chrom2 = infile.gettid("chr" + chrom2_info)

        if chrom == -1 or chrom2 == -1:
            logging.warning(f"Chromosome from record in --sites not found in input file header CHROM={r.chrom}, POS={r.start}, CHROM2={chrom2_info}, END={r.stop}")

        if isinstance(chrom, str):
            raise ValueError(f"Could not find {chrom} in bam file header")
        if isinstance(chrom2, str):
            chrom2 = chrom

        start = r.start  # note pysam converts to zero-based index like bam
        stop = r.stop
        if chrom == chrom2 and stop < start:
            start = r.stop - 1
            stop = r.start - 1

        if chrom not in recs:
            recs[chrom] = []

        if "SVLEN" in r.info:
            s = r.info["SVLEN"]
            if isinstance(s, tuple):
                svlen = abs(int(s[0]))
            else:
                svlen = int(r.info["SVLEN"])
        else:
            svlen = -1

        # is_dysgu = False
        # if "SVMETHOD" in r.info and "DYSGU" in r.info["SVMETHOD"]:
        #     is_dysgu = True
        #     seqs = parse_variant_seqs_dysgu(r, svt, r.chrom, start, chrom2_info, stop, paired_end, ref_genome)

        if parse_probs == "True":
            if "MeanPROB" in r.info:
                pval = r.info["MeanPROB"]
            elif len(r.samples) == 1:
                pval = r.samples[0].get("PROB")
                if pval is None:
                    pval = default_prob
                if pval == -1:
                    pval = default_prob
            else:
                pval = default_prob
        else:
            pval = default_prob

        if isinstance(pval, str) or pval < 0 or pval > 1:
            raise ValueError("PROB or MeanPROB must be in a float in range 0-1, error at CHROM {}, POS {}".format(chrom, start))

        recs[chrom].append(Site(chrom, start, chrom2, stop, svt, idx, r.id, svlen, pval))
    if not_parsed:
        logging.warning("Some records had incompatible SVTYPE: {}".format(
            str({k: v for k, v in Counter(not_parsed).items()})))
    return {k: deque(sorted(v, key=lambda x: x[1])) for k, v in recs.items()}


def append_uncalled(df, site_adder, infile, parse_probs):
    if site_adder is None:
        raise ValueError("Sites was None type")
    # add 0/0 genotype to uncalled variants in --sites
    found_sites = set([i for i in df["site_info"] if i])
    uncalled = []
    for chrom, sites in site_adder.sites.items():
        for k in sites:
            if k not in found_sites:
                r = {k: "." for k in df.columns}
                r["chrA"] = infile.get_reference_name(k.chrom)
                r["posA"] = k.start
                r["chrB"] = infile.get_reference_name(k.chrom2)
                r["posB"] = k.end
                r["svtype"] = k.svtype
                r["site_id"] = k.id
                r["svlen"] = k.svlen

                r["GT"] = "0/0"

                keys = ['GQ', 'DP', 'DN', 'DApri', 'DAsupp', 'NMpri', 'NMsupp', 'NMbase', 'MAPQpri', 'MAPQsupp', 'NP',
                        'maxASsupp', 'su', 'spanning', 'pe', 'supp', 'sc', 'bnd', 'sqc', 'scw', 'clip_qual_ratio',
                        'block_edge', 'raw_reads_10kb', 'mcov', 'linked', 'neigh', 'neigh10kb', 'ref_bases', "plus", "minus",
                        "strand_binom_t",
                        'n_gaps', "n_sa", "n_xa", "n_unmapped_mates", "double_clips", "remap_score", "remap_ed",
                        "bad_clip_count", "fcc", "n_small_tlen", "ras", 'fas', "inner_cn", "outer_cn", "compress",
                        "ref_rep", "jitter", 'prob']
                for fmt in keys:
                    r[fmt] = 0
                if parse_probs:
                    r["prob"] = k.prob
                r["gc"] = -1

                uncalled.append(r)

    df = pd.concat([df, pd.DataFrame.from_records(uncalled)])

    return df

