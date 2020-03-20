#cython: language_level=3

from __future__ import absolute_import
import multiprocessing
import sys
import pkg_resources
import click
import time
import datetime
import pickle
import numpy as np

from . import data_io, pairing, io_funcs


cdef void process_template(read_template):
    paired = io_funcs.sam_to_array(read_template)

    if paired:
        return

    res = pairing.process(read_template)

    if res:
        read_template["passed"] = True
        io_funcs.add_scores(read_template, *res)
        io_funcs.choose_supplementary(read_template)
        io_funcs.score_alignments(read_template, read_template["ri"], read_template['rows'], read_template['data'])


def load_mq_model(pth):
    if pth:
        click.echo("Loading MapQ recalibrator {}".format(pth), err=True)
        return pickle.load(open(pth, "rb"))
    else:
        click.echo("No MapQ recalibration", err=True)
        return None


def phred_from_model(p):
    if p == 1:
        return 30
    if p < 0.5:
        return 0
    # return 40
    P = 1 - p
    v = int(round(-10 * np.log10(P)))
    return v if v <= 30 else 30


def predict_mapq(xtest, model):
    return list(map(phred_from_model, model.predict_proba(xtest)[:, 1]))


cdef write_records(sam, mq_model, outsam):
    # Model features are "AS", "DA", "DN", "DP", "NP", "PS", "XS", "kind_key", "mapq", "DS
    if not mq_model:
        for name, record in sam:
            outsam.write(data_io.sam_to_str(name, record))

    else:
        map_qs = []
        for name, alns in sam:

            for a in alns:
                mapq = a[3]

                t = {i[:2]: i[5:] for i in a[10:]}
                kind_key = 0 if int(a[0]) & 2048 else 1
                f = [t["AS"], t["DA"], t["DN"], t["DP"], t["NP"], t["PS"], t["XS"] if "XS" in t else 0, kind_key, mapq, t["DS"]]
                map_qs.append(f)

        mqs = iter(predict_mapq(np.array(map_qs).astype(float), mq_model))

        for name, alns in sam:
            for i in range(len(alns)):
                flag = int(alns[i][0])
                # cigar = alns[i][4]
                mq = next(mqs)

                if flag & 2048:
                    if flag & 4:
                        alns[i][3] = "0"
                    else:
                        alns[i][3] = str(mq)

            outsam.write(data_io.sam_to_str(name, alns))


cpdef list job(data_tuple):

    temp = data_io.make_template(*data_tuple)

    process_template(temp)
    sam_temp = []
    if temp['passed']:
        sam = data_io.to_output(temp)
        if sam:
            sam_temp.append((temp["name"], sam))
    return sam_temp


def process_reads(args):
    t0 = time.time()

    click.echo("dysgu choose reading data from {}".format(args["sam"]), err=True)
    insert_std = args["template_size"].split(",")
    args["insert_median"] = float(insert_std[0])
    args["insert_stdev"] = float(insert_std[1])

    if not args["include"]:
        args["bias"] = 1.0
    else:
        click.echo("Elevating alignments in --include with --bias {}".format(args["bias"]), err=True)

    if args["output"] in {"-", "stdout"} or args["output"] is None:
        click.echo("Writing alignments to stdout", err=True)
        outsam = sys.stdout
    else:
        click.echo("Writing alignments to {}".format(args["output"]), err=True)
        outsam = open(args["output"], "w")

    if (args["fq1"] or args["fq2"]) and args["procs"] > 1:
        raise ValueError("Cant use procs > 1 with fq input")

    map_q_recal_model = load_mq_model(args["mq"])

    count = 0

    click.echo("dysgu {} process".format(args["procs"]), err=True)

    version = pkg_resources.require("dysgu")[0].version
    itr = data_io.iterate_mappings(args, version)

    isize = (args["insert_median"], args["insert_stdev"])
    match_score = args["match_score"]
    pairing_params = (args["max_insertion"], args["min_aln"], args["max_overlap"], args["ins_cost"],
                      args["ol_cost"], args["inter_cost"], args["u"])
    paired_end = int(args["paired"] == "True")
    bias = args["bias"]
    replace_hard = int(args["replace_hardclips"] == "True")

    max_d = args["insert_median"] + 4*args["insert_stdev"]  # Separation distance threshold to call a pair discordant

    if args["procs"] != 1:

        header_string = next(itr)
        outsam.write(header_string)

        temp = []
        for rows, last_seen_chrom, fq in itr:
            temp.append((rows, max_d, last_seen_chrom, fq, pairing_params, paired_end, isize,
                                         match_score, bias, replace_hard))

            if len(temp) > 10000:
                with multiprocessing.Pool(args["procs"]) as p:
                    res = p.map(job, temp)

                res = [item for sublist in res for item in sublist if item]  # Remove [] and flatten list
                write_records(res, map_q_recal_model, outsam)
                temp = []

        if len(temp) > 0:

            with multiprocessing.Pool(args["procs"]) as p:
                res = p.map(job, temp)
            res = [item for sublist in res for item in sublist if item]
            write_records(res, map_q_recal_model, outsam)

    # Use single process for debugging
    else:

        header_string = next(itr)
        outsam.write(header_string)

        sam_temp = []
        for rows, last_seen_chrom, fq in itr:

            count += 1

            # rows, max_d, last_seen_chrom, fq
            temp = data_io.make_template(rows, max_d, last_seen_chrom, fq, pairing_params, paired_end, isize,
                                         match_score, bias, replace_hard)

            process_template(temp)

            if temp['passed']:
                sam = data_io.to_output(temp)
                if sam:
                    sam_temp.append((temp["name"], sam))

                if len(sam_temp) > 50000:
                    write_records(sam_temp, map_q_recal_model, outsam)
                    sam_temp = []

        if len(sam_temp) > 0:
            write_records(sam_temp, map_q_recal_model, outsam)

    if args["output"] != "-" or args["output"] is not None:
        outsam.close()

    click.echo("dysgu choose {} completed in {} h:m:s".format(args["sam"],
                                                            str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)
