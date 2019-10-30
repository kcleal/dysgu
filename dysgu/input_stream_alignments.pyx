#cython: language_level=3

from __future__ import absolute_import
import multiprocessing
import sys
import pkg_resources
from threading import Thread
import click
from dysgu import data_io, pairing, io_funcs
import time
import datetime
import pickle
import numpy as np


def process_template(read_template):
    # click.echo(read_template, err=True)
    # quit()
    # if read_template["name"] == "HISEQ1:11:H8GV6ADXX:1:2113:3986:26065":
    #     click.echo(read_template, err=True)
    #     quit()

    paired = io_funcs.sam_to_array(read_template)

    if paired:
        return

    res = pairing.process(read_template)
    # click.echo(res, err=True)
    if res:
        read_template["passed"] = True
        io_funcs.add_scores(read_template, *res)
        io_funcs.choose_supplementary(read_template)
        io_funcs.score_alignments(read_template, read_template["ri"], read_template['rows'], read_template['data'])


def worker(queue, out_queue):
    while True:
        job = queue.get(True)

        if job == "Done":
            queue.task_done()
            out_queue.put("Done")
            break

        else:
            big_string = ""
            for data_tuple in job:
                read_template = data_io.make_template(*data_tuple)
                process_template(read_template)
                if read_template['passed']:
                    outstring = data_io.to_output(read_template)
                    if outstring:
                        big_string += outstring
                    else:
                        pass  # No mappings

            if len(big_string) > 0:
                out_queue.put(big_string)
            # else:
            #     click.echo("WARNING: no output from job.", err=True)
        queue.task_done()


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


def write_records(sam, mq_model, outsam):
    # Model features are "AS", "DA", "DN", "DP", "NP", "PS", "XS", "kind_key", "mapq", "DS
    if not mq_model:
        for name, record in sam:
            outsam.write(data_io.sam_to_str(name, record))

    else:  # Todo multiprocessing here, cython for array preparation etc
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


def process_reads(args):
    t0 = time.time()

    click.echo("dysgu choose reading data from {}".format(args["sam"]), err=True)

    if not args["include"]:
        args["bias"] = 1.0
    else:
        click.echo("Elevating alignments in --include with --bias {}".format(args["bias"]), err=True)

    if args["output"] == "-" or args["output"] is None:
        click.echo("Writing alignments to stdout", err=True)
        outsam = sys.stdout
    else:
        click.echo("Writing alignments to {}".format(args["output"]), err=True)
        outsam = open(args["output"], "w")

    map_q_recal_model = load_mq_model(args["mq"])

    count = 0

    click.echo("fnfi {} process".format(args["procs"]), err=True)

    version = pkg_resources.require("dysgu")[0].version
    itr = data_io.iterate_mappings(args, version)

    # Use multiprocessing:
    # https://stackoverflow.com/questions/17241663/filling-a-queue-and-managing-multiprocessing-in-python
    if args["procs"] != 1:

        cpus = args["procs"] if args["procs"] != 0 else multiprocessing.cpu_count()
        click.echo("fnfi align runnning {} cpus".format(cpus), err=True)

        # Todo joinable queue is not efficient, other options?
        the_queue = multiprocessing.JoinableQueue(maxsize=100)  #cpus+2)
        out_queue = multiprocessing.Queue()

        the_pool = multiprocessing.Pool(args["procs"] if args["procs"] != 0 else multiprocessing.cpu_count(),
                                        worker, (the_queue, out_queue,))

        def writer_thread(q, outsam):
            while True:
                aln = q.get()
                if aln == "Done":
                    break
                elif aln == "Job failed":
                    click.echo("job failed", err=True)
                elif len(aln) > 1:
                    outsam.write(aln)

            click.echo("Writing done", err=True)

        writer = Thread(target=writer_thread, args=(out_queue, outsam, ))
        writer.setDaemon(True)
        writer.start()

        job = []
        header_string = next(itr)
        out_queue.put(header_string)

        for data_tuple in itr:
            count += 1

            job.append(data_tuple)
            if len(job) > 500:
                the_queue.put(job)
                job = []
        if len(job) > 0:
            the_queue.put(job)

        the_queue.join()  # Wait for jobs to finish
        the_queue.put("Done")  # Send message to stop workers
        writer.join()  # Wait for writer to closing

    # Use single process for debugging
    else:

        header_string = next(itr)
        outsam.write(header_string)

        sam_temp = []
        for data_tuple in itr:

            count += 1
            temp = data_io.make_template(*data_tuple)

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

    # if args["output"] != "-" or args["output"] is not None:
    #     outsam.close()

    click.echo("dysgu choose {} completed in {} h:m:s".format(args["sam"],
                                                            str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)


if __name__ == "__main__":

    t = {'isize': (210.0, 175.0), 'max_d': 910.0, 'match_score': 1.0, 'pairing_params': (150.0, 17.0, 150.0, 0.1, 1.0, 2.0, 9.0), 'paired_end': 1, 'inputdata': [(['chr16:53948053-53948078+chr16:46380752-46380852_chr16:46380794-46380919', '99', 'chr16', '46380753', "60\t25S100M\t=\t46380795\t167\tGAACTCGCCATTGTAAAGTAAGAGAATTCCATTCAATTCCATTTGATGATGTTTCTCTTGTATTCCATTGGATAATTCCTTTCAGTTCCCTACGATGATTATTCTTTCGAGTCAATTCGATGATC\tGGEGGD;GFD5F?;2GGGGGG-FFGFEG9GFFC;FGGGFGC;GGGG5GDGDGGFFGFF2FG@FGGFFGA;GGGEGFGG-BGGF@GFGGGA@DFGGFGFFG@8D3FDGGDFGGD@1?0GGGGFA'$\tNM:i:1\tMD:Z:99T0\tAS:i:99\tXS:i:79\n"], 0), (['chr16:53948053-53948078+chr16:46380752-46380852_chr16:46380794-46380919', '147', 'chr16', '46380795', '60\t125M\t=\t46380753\t-167\tTTGGATAATTCCTTTCAGTTCCCTACGCTGATTATTCTTTCGAGTCAATTCGATGATTCTATTCCATTCCCTTCGATGATGATTCCATTTCACTCCATTTGATGATTCCATTCGACTCAATTTGG\tD=;4=E5GD015DGGBDGFFBFF;FEA9FGFGGGFGGCGF>E>FGFAD?FFCFGEDGGBGGGGF@GG<GGGGGGGGFFFFEEECCBGGG7=GGF2?FFGBAFGGGGGGFGGGGGGGGFGGGGGGF\tNM:i:1\tMD:Z:27A97\tAS:i:120\tXS:i:100\n'], 0)], 'bias': 1.0, 'read1_length': 0, 'read2_length': 0, 'score_mat': {}, 'passed': 0, 'name': 'chr16:53948053-53948078+chr16:46380752-46380852_chr16:46380794-46380919', 'last_seen_chrom': 'chr7', 'inputfq': (None, None), 'read1_seq': 0, 'read2_seq': 0, 'read1_q': 0, 'read2_q': 0, 'read1_reverse': 0, 'read2_reverse': 0, 'replace_hard': 0, 'fq_read1_seq': 0, 'fq_read2_seq': 0, 'fq_read1_q': 0, 'fq_read2_q': 0}

    # print(t.keys())

    process_template(t)
    print(t.keys())
    print(t["passed"])
    print(t["outstr"])
    print(t["read1_seq"])
    print(t["read2_seq"])
    print(t["read1_length"])
    print(t["read2_length"])
    print(t["score_mat"])
    # print(t["rows"])

    sam = data_io.to_output(t)

    for item in sam:
        print(item)

    print()

    for item in t["inputdata"]:
        print(item)