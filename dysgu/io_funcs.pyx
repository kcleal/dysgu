#!python
#cython: language_level=2, boundscheck=False, wraparound=False
#distutils: language=c++


import numpy as np
cimport numpy as np
cimport cython
import click
from collections import defaultdict
import ncls


DTYPE = np.float
ctypedef np.float_t DTYPE_t


from libc.stdlib cimport malloc
import re

cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0',
                       '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
                       '\0', '\0', '\0', '\0',  'a',  'a' ]

np.random.seed(0)


def reverse_complement(str seq, int seq_len):
    """https://bioinformatics.stackexchange.com/questions/3583/\
    what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho/3595#3595"""

    cdef char *seq_dest = <char *>malloc(seq_len + 1)
    seq_dest[seq_len] = '\0'

    cdef bytes py_bytes = seq.encode('UTF-8')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('UTF-8')


def get_start_end(str cigar):
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    cdef int end = 0
    cdef int start = 0
    cdef int i

    for i in range(0, len(c)-1, 2):
        if i == 0 and (c[i+1] == "S" or c[i+1] == "H"):
            start += int(c[i])
            end += int(c[i])
        elif c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += int(c[i])
    return start, end


def get_align_end_offset(str cigar):
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    cdef int end = 0
    cdef int i
    for i in range(0, len(c)-1, 2):
        if c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += int(c[i])
    return end


def check_for_good_pairing(r1, r2, name, max_d):
    # Check to see if this pair of alignments needs pairing
    cdef int aflag = int(r1[0])
    cdef int bflag = int(r2[0])

    # Check if rnext and pnext set properly
    cdef int p1, p2
    cdef int proper_pair = 0
    cdef str dn = '250.0'
    if r1[2] == r2[6] and r2[2] == r1[6]:  # Might not generalize for mappers other than bwa mem

        p1 = int(r1[2])
        p2 = int(r2[2])

        if (aflag & 16 and not bflag & 16) or (not aflag & 16 and bflag & 16):  # Not on same strand

            if abs(p1 - p2) < max_d:
                # Check for FR or RF orientation
                if (p1 < p2 and (not aflag & 16) and (bflag & 16)) or (p2 <= p1 and (not bflag & 16) and (aflag & 16)):
                    proper_pair = 1

        if proper_pair == 1:
            dn = '0.0'
        # Note DP and PS are clipped at 250
        tags = ["SP:Z:0", "DA:i:100", "DP:Z:250.0", "DN:Z:" + dn, "PS:Z:250.0", "NP:Z:1.0", "DS:i:0"]
        r1 += tags
        r2 += tags

        return [r1, r2]

    return 0


def sam_to_array(template):
    # Expect read1 and read2 alignments to be concatenated, not mixed together
    data, overlaps = list(zip(*template["inputdata"]))
    template["inputdata"] = [[i[1], i[2], i[3]] + i[4].strip().split("\t") for i in data]

    # If only one alignment for read1 and read2, no need to try pairing, just send sam to output
    if len(data) == 2:
        pair_str = check_for_good_pairing(template["inputdata"][0], template["inputdata"][1], template["name"], template["max_d"])

        if pair_str:
            template["passed"] = 1
            template["outstr"] = pair_str
            return 1

    # [chrom, pos, query_start, query_end, aln_score, row_index, strand, read, num_mis-matches, original_aln_score]
    cdef np.ndarray[np.float_t, ndim=2] arr = np.zeros((len(data), 10))  # , dtype=np.float

    chrom_ids = {}

    if template["inputdata"][0][1] == "*":
        template["inputdata"][0][1] = template["last_seen_chrom"]

    cdef int cc = 0
    cdef int idx, pos, flag, seq_len, query_start, query_end, start_temp  # new_start, new_end,
    cdef str chromname, cigar, k, t, v
    cdef float bias = template["bias"]

    cdef int read1_strand_set, read2_strand_set, current_l
    cdef int first_read2_index = len(template["inputdata"])
    read1_set = 0  # Occationally multiple primaries, take the longest
    read2_set = 0

    for idx in range(len(template["inputdata"])):

        l = template["inputdata"][idx]

        if l[1] != "*":
            chromname = l[1]
            template["last_seen_chrom"] = chromname
        else:
            l[1] = template["last_seen_chrom"]
            chromname = l[1]
        if chromname not in chrom_ids:
            chrom_ids[chromname] = cc
            cc += 1
        arr[idx, 0] = chrom_ids[chromname]  # l.rname  # chrom name

        pos = int(l[2])  # Add hard clips
        arr[idx, 1] = pos  # l.pos
        arr[idx, 5] = idx
        flag = int(l[0])
        arr[idx, 6] = -1 if flag & 16 else 1  # Flag

        # Sometimes both first and second are not flagged. Assume first
        if not flag & 64 and not flag & 128:
            arr[idx, 7] = 1
        else:
            if flag & 64:
                arr[idx, 7] = 1
            else:
                arr[idx, 7] = 2
                if idx < first_read2_index:
                    first_read2_index = idx

        tags = [i.split(":") for i in l[11:]]
        seq_len = len(l[8])

        for k, t, v in tags:
            if k == "NM":
                arr[idx, 8] = float(v)
            elif k == "AS":
                arr[idx, 9] = float(v)  # Keep copy of original alignment score
                if overlaps[idx]:
                    arr[idx, 4] = float(v) * bias
                else:
                    arr[idx, 4] = float(v)

        current_l = len(l[9])

        if template["paired_end"]:
            if flag & 64 and read1_set < current_l and len(l[8]) > 1:  # First in pair
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]
                template["read1_length"] = seq_len
                read1_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template["read1_reverse"] = 1 if flag & 16 else 0

            elif flag & 128 and read2_set < current_l and len(l[8]) > 1:  # First in pair
                template["read2_seq"] = l[8]
                template["read2_q"] = l[9]
                template["read2_length"] = seq_len
                read2_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template["read2_reverse"] = 1 if flag & 16 else 0

        else:
            if template["read1_seq"] == 0 and not (flag & 256) and (len(l[8]) > 1) and read1_set < current_l:
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]
                template["read1_length"] = len(l[8])
                read1_set = current_l

        cigar = l[4]
        if not cigar:
            query_start = 0  # Unmapped read, no cigar
            query_end = 0
        else:
            query_start, query_end = get_start_end(cigar)

            # If current alignment it not primary, and on different strand from primary, count from other end
            if template["paired_end"]:
                # if flag & 256 or flag & 2048:  # Todo is this right? Are primary reads counted from end too?
                    if flag & 64 and template["read1_reverse"] != bool(flag & 16):
                        start_temp = template["read1_length"] - query_end
                        query_end = start_temp + query_end - query_start
                        query_start = start_temp
                            # Different strand to primary, count from end
                            # new_end = template["read1_length"] - query_start
                            # new_start = template["read1_length"] - query_end
                            # query_start = new_start
                            # query_end = new_end

                    elif flag & 128 and (template["read2_reverse"] != bool(flag & 16)):
                        # new_end = template["read2_length"] - query_start
                        # new_start = template["read2_length"] - query_end
                        # query_start = new_start
                        # query_end = new_end
                        start_temp = template["read2_length"] - query_end
                        query_end = start_temp + query_end - query_start
                        query_start = start_temp

            if not template["paired_end"]:
                # click.echo(("qstart", query_start, query_end, template["read1_length"]), err=True)
                if flag & 16:  # Single end Reverse strand, count from end
                    start_temp = template["read1_length"] - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp



        # if template["name"] == "1-WT-parent_length_367_mean_cov_5.17711171662":
        #     click.echo((int(query_start), int(query_end), flag, template["paired_end"]), err=True)
        arr[idx, 2] = query_start
        arr[idx, 3] = query_end  # query_start + query_end


    # Save any input fastq information
    fq1, fq2 = template["inputfq"]

    if fq1:
        template["fq_read1_seq"] = fq1[1]
        template["fq_read1_q"] = fq1[2]
        template["read1_length"] = len(fq1[1])
    if fq2:
        template["fq_read2_seq"] = fq2[1]
        template["fq_read2_q"] = fq2[2]
        template["read2_length"] = len(fq2[1])

    cdef int j
    if template["paired_end"]:  # Increment the contig position of read 2

        for j in range(first_read2_index, len(arr)):
            if arr[j, 7] == 2:  # Second in pair
                arr[j, 2] += template['read1_length']
                arr[j, 3] += template['read1_length']

            if arr[j, 3] > template["read1_length"] + template["read2_length"]:
                click.echo((template["read1_length"], template["read2_length"]), err=True)
                click.echo(arr[j, 3], err=True)
                raise ValueError
    #np.random.shuffle(arr)  # Randomize order of chromosomes
    template["first_read2_index"] = first_read2_index
    # if template["name"] == "1-WT-parent_length_367_mean_cov_5.17711171662":
    #     click.echo(template["name"], err=True)
    #     for roww in arr:
    #         click.echo(list(roww.astype(int)), err=True)
    template['data'] = np.array(sorted(arr, key=lambda x: (x[2], -x[4]))).astype(float)
    template['chrom_ids'] = chrom_ids
    # quit()
    del template["inputfq"]

    return 0


def choose_supplementary(template):

    cdef int j = 0
    template['ri'] = dict(zip(template['data'][:, 5], range(len(template['data']))))  # Map of row_index and array index

    cdef np.ndarray[long, ndim=1] actual_rows = np.array([template['ri'][j] for j in template['rows']]).astype(int)
    cdef np.ndarray[double, ndim=2] d = template['data']  #[actual_rows, :]  # np.float_t is double

    cdef double read1_max = 0
    cdef double read2_max = 0
    cdef int i = 0

    for j in range(len(actual_rows)):
        i = actual_rows[j]
        if d[i, 7] == 1 and d[i, 9] > read1_max:  # Use original alignment score, not biased
            read1_max = d[i, 9]

        elif d[i, 7] == 2 and d[i, 9] > read2_max:
            read2_max = d[i, 9]

    ids_to_name = {v: k for k, v in template["chrom_ids"].items()}

    locs = []
    cdef double m = 0
    for j in range(len(actual_rows)):
        i = actual_rows[j]

        loc = "{}-{}-{}-{}".format(ids_to_name[int(d[i, 0])], int(d[i, 1]), int(d[i, 6]), int(d[i, 7]) )
        locs.append(loc)

        if loc not in template['score_mat']:
                template['score_mat'][loc] = []
        # Values are popped when setting supplementary; prevents bug where read contains two identical aligns
        if d[i, 7] == 1:
            m = read1_max
        else:
            m = read2_max

        if d[i, 9] == m:  # Primary, next best s
            template['score_mat'][loc] += [True, 0]
        else:
            template['score_mat'][loc] += [False, 0]
    template['locs'] = locs


def score_alignments(template, ri, np.ndarray[np.int64_t, ndim=1]  template_rows, np.ndarray[DTYPE_t, ndim=2] template_data):
    # Scans all alignments for each query, slow for long reads but ok for short read data

    all_xs = []
    cdef int i, actual_row, item, idx
    cdef float xs = -1
    cdef float size = 0
    cdef float qstart = 0
    cdef float qend = 0
    cdef float isize = 0
    cdef float istart = 0
    cdef float iend = 0
    cdef float iscore = 0
    cdef float ol = 0
    cdef float ori_aln_score = 0
    cdef int twoi = template["first_read2_index"]  # Use to skip read1/read2 alignments

    idx = 0
    for item in template_rows:
        actual_row = ri[item]
        qstart, qend, readn, ori_aln_score = template_data[actual_row, [2, 3, 7, 9]]
        size = qend - qstart

        if template["paired_end"]:
            if readn == 2:
                rr = range(twoi, len(template_data))
            else:
                rr = range(0, twoi)
        else:
            rr = range(len(template_data))

        for i in rr:
            if i == actual_row:
                continue
            istart, iend, iscore = template_data[i, [2, 3, 4]]  # Note use biased align score, otherwise gets confusing
            isize = (iend - istart) + 1e-6
            # Check for overlap
            ol = max(0, min(qend, iend) - max(qstart, istart))
            if ol and (ol / size > .85) and (ol / isize > .85):  # check for 85% reciprocal overlap with current
                if iscore > xs:
                    xs = iscore

        if xs == -1:
            xs = ori_aln_score
        template["score_mat"][template["locs"][idx]][1] = xs
        idx += 1


def add_scores(template, np.ndarray[np.float_t, ndim=1] rows, float path_score, float second_best, float dis_to_normal, int norm_pairings):

    # The rows correspond to the indexes of the input array, not the original ordering of the data
    template['rows'] = rows.astype(int)  # list(map(int, rows))

    # To get the original rows use, col 5 is the row index: mapped to actual index
    template['score_mat']["dis_to_next_path"] = path_score - second_best
    template['score_mat']["dis_to_normal"] = dis_to_normal
    template['score_mat']["path_score"] = path_score
    template['score_mat']['normal_pairings'] = norm_pairings


cdef list get_bed_regions(str bed):
    b = [tuple([int(j) if j.isdigit() else j for j in i.strip().split("\t")[:3]]) for i in open(bed, "r")
         if i[0] != "#" and len(i) > 0 and "\t" in i]
    if len(b) == 0:
        raise ValueError("Bed regions not formatted correctly")
    return b


cpdef dict overlap_regions(str bed, int_chroms=False, infile=None):
    if not bed:
        return {}
    regions = get_bed_regions(bed)
    chrom_interval_start = defaultdict(list)
    chrom_interval_end = defaultdict(list)
    for c, s, e in regions:
        if int_chroms:
            c = infile.gettid(c)
        chrom_interval_start[c].append(int(s))
        chrom_interval_end[c].append(int(e))

    regions = {k: ncls.NCLS(np.array(chrom_interval_start[k]),
                            np.array(chrom_interval_end[k]),
                            np.array(chrom_interval_start[k])) for k in chrom_interval_start}

    return regions


cpdef int intersecter_int_chrom(dict tree, int chrom, int start, int end):
    if not tree:
        return 0
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return 1
        else:
            return 0
    else:
        return 0


cpdef int intersecter_str_chrom(dict tree, str chrom, int start, int end):
    if not tree:
        return 0
    elif chrom in tree:
        if len(list(tree[chrom].find_overlap(start, end))) > 0:
            return 1
        else:
            return 0
    else:
        return 0