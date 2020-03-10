#!python
#cython: language_level=2, boundscheck=False, wraparound=False
#distutils: language=c++


import numpy as np
cimport numpy as np
cimport cython
import click
from collections import defaultdict
import ncls
import pkg_resources
import sortedcontainers
import pandas as pd


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


cpdef str reverse_complement(str seq, int seq_len):
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


cpdef tuple get_start_end(str cigar):
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


cpdef int get_align_end_offset(str cigar):
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    cdef int end = 0
    cdef int i
    for i in range(0, len(c)-1, 2):
        if c[i+1] not in "DHS":  # Don't count deletions, or soft/hard clips at right-hand side
            end += int(c[i])
    return end


cdef tuple check_for_good_pairing(template):

    # Check to see if this pair of alignments needs pairing
    r1 = template["inputdata"][0]
    r2 = template["inputdata"][1]

    cdef int aflag = int(r1[0])
    cdef int bflag = int(r2[0])

    if aflag & 4:
        template["read1_unmapped"] = 1
    if bflag & 4:
        template["read2_unmapped"] = 1

    # Check if rnext and pnext set properly
    cdef int p1, p2
    cdef int proper_pair = 1 if aflag & 2 else 0
    cdef str dn = '250.0'
    # if r1[2] == r2[6] and r2[2] == r1[6]:  # Might not generalize for mappers other than bwa mem

    if not proper_pair and not aflag & 12:  # Read unmapped/mate unmapped. Double check in case flag not set properly
        p1 = int(r1[2])  # Position
        p2 = int(r2[2])

        if (aflag & 16 and not bflag & 16) or (not aflag & 16 and bflag & 16):  # Not on same strand

            if abs(p1 - p2) < template["max_d"]:
                # Check for FR or RF orientation
                if (p1 < p2 and (not aflag & 16) and (bflag & 16)) or (p2 <= p1 and (not bflag & 16) and (aflag & 16)):
                    proper_pair = 1

    if proper_pair == 1:
        dn = '0.0'
    else:
        dn = '250.0'
    # Note DP and PS are clipped at 250
    # tags = ["SP:Z:0", "DA:i:100", "DP:Z:250.0", "DN:Z:" + dn, "PS:Z:250.0", "NP:Z:1.0", "DS:i:0"]

    if aflag & 4 and bflag & 4:  # Read unmapped, mate unmapped
        r1 += ["SP:Z:0", "DA:i:0", "DP:Z:0", "DN:Z:" + dn, "PS:Z:0", "NP:Z:0", "DS:i:0"]
        r2 += ["SP:Z:0", "DA:i:0", "DP:Z:0", "DN:Z:" + dn, "PS:Z:0", "NP:Z:0", "DS:i:0"]

    elif aflag & 4:
        r1 += ["SP:Z:0", "DA:i:0", "DP:Z:0", "DN:Z:" + dn, "PS:Z:0", "NP:Z:0", "DS:i:0"]
        r2 += ["SP:Z:0", "DA:i:100", "DP:Z:125", "DN:Z:" + dn, "PS:Z:125", "NP:Z:0", "DS:i:0"]

    elif bflag & 4:
        r1 += ["SP:Z:0", "DA:i:100", "DP:Z:125", "DN:Z:" + dn, "PS:Z:125", "NP:Z:0", "DS:i:0"]
        r2 += ["SP:Z:0", "DA:i:0", "DP:Z:0", "DN:Z:" + dn, "PS:Z:0", "NP:Z:0", "DS:i:0"]

    else:
        r1 += ["SP:Z:0", "DA:i:100", "DP:Z:250", "DN:Z:" + dn, "PS:Z:250", "NP:Z:1", "DS:i:0"]
        r2 += ["SP:Z:0", "DA:i:100", "DP:Z:250", "DN:Z:" + dn, "PS:Z:250", "NP:Z:1", "DS:i:0"]

    # r1 += tags
    # r2 += tags

    return r1, r2


def sort_func(row):
    return row[7], row[2], -row[4]


cpdef sam_to_array(template):
    # Expect read1 and read2 alignments to be concatenated, not mixed together
    data, overlaps = list(zip(*template["inputdata"]))
    template["inputdata"] = [[i[1], i[2], i[3]] + i[4].strip().split("\t") for i in data]

    # If only one alignment for read1 and read2, no need to try pairing, just send sam to output
    if template["paired_end"] and len(data) == 2:
        pair_str = check_for_good_pairing(template)# template["inputdata"][0], template["inputdata"][1], template["name"], template["max_d"])

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
    cdef int first_read2_index = len(template["inputdata"]) + 1
    read1_set = 0  # Occasionally multiple primaries, take the longest
    read2_set = 0

    srt_1 = []
    srt_2 = []

    for idx in range(len(template["inputdata"])):

        l = template["inputdata"][idx]

        flag = int(l[0])
        pos = int(l[2])  # Add hard clips

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
        arr[idx, 1] = pos  # l.pos
        arr[idx, 5] = idx
        arr[idx, 6] = -1 if flag & 16 else 1  # Flag

        if idx == 0 and flag & 4:
            template["read1_unmapped"] = 1

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
                    if flag & 4:
                        template["read2_unmapped"] = 1

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
                srt_2.append(-float(v))
        current_l = len(l[9])

        if template["paired_end"]:

            if flag & 64 and read1_set < current_l and len(l[8]) > 1:  # First in pair
                template["read1_seq"] = l[8]
                template["read1_q"] = l[9]
                template["read1_length"] = seq_len
                read1_set = current_l
                if not (flag & 256 or flag & 2048):  # Set primary read strand
                    template["read1_reverse"] = 1 if flag & 16 else 0

            elif flag & 128 and read2_set < current_l and len(l[8]) > 1:  # Second in pair

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
            query_start = 0  # Unmapped read? no cigar
            query_end = 0
            srt_1.append(query_start)

        else:
            query_start, query_end = get_start_end(cigar)
            srt_1.append(query_start)

            # If current alignment it not primary, and on different strand from primary, count from other end
            if template["paired_end"]:
                if flag & 64 and template["read1_reverse"] != bool(flag & 16):  # First in pair, read1_rev != read_rev
                    start_temp = template["read1_length"] - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp

                elif flag & 128 and (template["read2_reverse"] != bool(flag & 16)):  # Second in pair
                    start_temp = template["read2_length"] - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp

            if not template["paired_end"]:
                if flag & 16:  # Single end Reverse strand, count from end
                    start_temp = template["read1_length"] - query_end
                    query_end = start_temp + query_end - query_start
                    query_start = start_temp


        arr[idx, 2] = query_start
        arr[idx, 3] = query_end  # query_start + query_end

    if first_read2_index == len(arr) + 1:
        template["paired_end"] = 0

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
                raise ValueError

    template["first_read2_index"] = first_read2_index

    template['data'] = np.array(sorted(arr, key=sort_func))

    template['chrom_ids'] = chrom_ids

    del template["inputfq"]

    return 0


cpdef choose_supplementary(dict template):
    # Final alignments have been chosen, but need to decide which is supplementary
    cdef int j = 0
    template['ri'] = dict(zip(template['data'][:, 5], range(len(template['data']))))  # Map of row_index and array index
    cdef np.ndarray[long, ndim=1] actual_rows
    try:
        actual_rows = np.array([template['ri'][j] for j in template['rows']]).astype(int)
    except:
        click.echo(template["ri"], err=True)
        click.echo(template["rows"], err=True)
        click.echo(template["inputdata"], err=True)
        click.echo((template["read1_unmapped"], template["read2_unmapped"]), err=True)
        quit()
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


cpdef void score_alignments(dict template, ri, np.ndarray[np.int64_t, ndim=1]  template_rows, np.ndarray[DTYPE_t, ndim=2] template_data):
    # Scans all alignments for each query, slow for long reads but ok for short read data
    # Used for DN, similar to XS
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


cpdef list col_names():
    return ["chrA", "posA", "chrB", "posB", "sample", "id", "kind", "svtype", "join_type", "cipos95A", "cipos95B",
         "DP", "DN", "DApri", "DAsupp",  "NMpri", "NMsupp", "MAPQpri", "MAPQsupp", "NP",
          "maxASsupp",  "su", "pe", "supp", "sc", "block_edge",
         "raw_reads_10kb",
          "linked", "contigA", "contigB",  "gc", "neigh", "rep", "rep_sc", "ref_bases", "svlen", "plus", "minus",
            ]  # "Prob"


def make_main_record(r, version, index, format_f, df_rows, add_kind):

    # Pick best row (best support, or highest prob if available
    if len(format_f) > 1:

        best = sorted([(int(v["su"]), k) for k, v in df_rows.items()], reverse=True)[0][1]
        r = df_rows[best]
        gc = r["gc"]
        rep = r["rep"]
        repsc = r["rep_sc"]
        su, pe, sr, sc = 0, 0, 0, 0
        # probs = []
        for row in df_rows.values():
            pe += row["pe"]
            sr += row["supp"]
            sc += row["sc"]
            su += (row["pe"] + row["supp"])
            # probs.append(row["Prob"])
        # probs = round(np.median(probs), 3)

    else:
        pe = r["pe"]
        sr = r["supp"]
        sc = r["sc"]
        su = (r["pe"] + r["supp"])
        # probs = r["Prob"]
        gc = r["gc"]
        rep = r["rep"]
        repsc = r["rep_sc"]

    samp = r["sample"]

    if r["chrA"] == r["chrB"] and r["posA"] > r["posB"]:
        chrA, posA, cipos95A, contig2 = r["chrA"], r["posA"], r["cipos95A"], r["contigB"]
        r["chrA"] = r["chrB"]
        r["posA"] = r["posB"]
        r["cipos95A"] = r["cipos95B"]
        r["chrB"] = chrA
        r["posB"] = posA
        r["cipos95B"] = cipos95A
        r["contigB"] = r["contigA"]
        r["contigA"] = contig2

    info_extras = []
    if r["chrA"] == r["chrB"]:
        svlen = abs(r["posA"] - r["posB"])
        info_extras.append(f"SVLEN={svlen}")

    if r["contigA"]:
        info_extras.append(f"CONTIGA={r['contigA']}")
    if r["contigB"]:
        info_extras.append(f"CONTIGB={r['contigB']}")

    if add_kind:
        info_extras += [f"KIND={r['kind']}"]

    info_extras += [f"GC={gc}",
                    f"REP={'%.3f' % rep}",
                    f"REPSC={'%.3f' % repsc}",
                    f"SU={su}",
                    f"PE={pe}",
                    f"SR={sr}",
                    f"SC={sc}",]
                    #f"MPROB={probs}"]

    rec = [r["chrA"], r["posA"], index, ".", f"<{r['svtype']}>", ".", ".",
           # INFO line
           ";".join([f"SVMETHOD=DYSGUv{version}",
                   f"SVTYPE={r['svtype']}",
                   f"END={r['posB']}",
                   f"CHR2={r['chrB']}",
                   f"CT={r['join_type']}",
                   f"CIPOS95={r['cipos95A']}",
                   f"CIEND95={r['cipos95B']}",

                   ] + info_extras),
           "GT:DP:DN:DAP:DAS:NMP:NMS:MAPQP:MAPQS:NP:MAS:SU:PE:SR:SC:BE:COV:LNK:NEIGH:RB:PS:MS"  # :PROB
           ]
    # FORMAT line(s)
    for item in format_f.values():
        rec.append(":".join(map(str, item)))

    return rec


def gen_format_fields(r, df, names):

    if len(names) == 1:
        return {0: (["./.", r['DP'], r['DN'], r['DApri'], r['DAsupp'], r['NMpri'], r['NMsupp'], r['MAPQpri'],
                                  r['MAPQsupp'], r['NP'], r['maxASsupp'], r['pe'] + r['supp'], r['pe'], r['supp'],
                                  r['sc'], r['block_edge'], r['raw_reads_10kb'], r['linked'], r['neigh'],
                                  r['ref_bases'], r["plus"], r["minus"]])}, {}

    cols = {}
    if "partners" in r:
        if not isinstance(r["partners"], list):
            if len(r["partners"]) == 0 or pd.isna(r["partners"]):
                r["partners"] = []
            else:
                r["partners"] = [int(i.split(",")[1]) for i in r["partners"].split("|")]

        for idx in r["partners"]:
            r2 = df.loc[idx]  # iloc
            cols[r2["table_name"]] = r2

    if "table_name" in r:
        cols[r["table_name"]] = r

    format_fields = sortedcontainers.SortedDict()

    for name in names:

        if name in cols:
            format_fields[name] = (["./.", r['DP'], r['DN'], r['DApri'], r['DAsupp'], r['NMpri'], r['NMsupp'], r['MAPQpri'],
                                  r['MAPQsupp'], r['NP'], r['maxASsupp'], r['pe'] + r['supp'], r['pe'], r['supp'],
                                  r['sc'], r['block_edge'], r['raw_reads_10kb'], r['linked'], r['neigh'],
                                  r['ref_bases'], r["plus"], r["minus"]])  # r['Prob']
        else:
            format_fields[name] = [0] * 20

    return format_fields, cols



def to_vcf(df, args, names, outfile, show_names=True,  contig_names=""):

    HEADER = """##fileformat=VCFv4.2
##source=DYSGU
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=CIPOS95,Number=1,Type=Integer,Description="Confidence interval size (95%) around POS for imprecise variants">
##INFO=<ID=CIEND95,Number=1,Type=Integer,Description="Confidence interval size (95%) around END for imprecise variants">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=KIND,Number=1,Type=String,Description="Kind of join with respect to input regions">
##INFO=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of supplementary reads supporting the variant across all samples">
##INFO=<ID=SC,Number=1,Type=Integer,Description="Number of soft-clip reads supporting the variant across all samples">
##INFO=<ID=CONTIGA,Number=1,Type=String,Description="Contig from CHROM POS">
##INFO=<ID=CONTIGB,Number=1,Type=String,Description="Contig from CHR2 END">
##INFO=<ID=GC,Number=1,Type=Float,Description="GC% of assembled contigs">
##INFO=<ID=REP,Number=1,Type=Float,Description="Repeat score for contigs aligned bases">
##INFO=<ID=REPSC,Number=1,Type=Float,Description="Repeat score for contigs soft-clipped bases">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean distance-to-pair metric supporting the variant">
##FORMAT=<ID=DN,Number=1,Type=Float,Description="Mean distance-to-normal metric supporting the variant">
##FORMAT=<ID=DAP,Number=1,Type=Float,Description="Mean distance-to-alignment metric for primary alignments">
##FORMAT=<ID=DAS,Number=1,Type=Float,Description="Mean distance-to-alignment metric for supplementary alignments">
##FORMAT=<ID=NMP,Number=1,Type=Float,Description="Mean edit distance for primary alignments supporting the variant">
##FORMAT=<ID=NMS,Number=1,Type=Float,Description="Mean edit distance for supplementary alignments supporting the variant">
##FORMAT=<ID=MAPQP,Number=1,Type=Float,Description="Mean MAPQ for primary reads supporting the variant">
##FORMAT=<ID=MAPQS,Number=1,Type=Float,Description="Mean MAPQ for supplementary reads supporting the variant">
##FORMAT=<ID=NP,Number=1,Type=Integer,Description="Number of alignments in normal-pair orientation supporting the variant">
##FORMAT=<ID=MAS,Number=1,Type=Integer,Description="Maximum alignment score of supplementary reads supporting the variant">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supplementary alignments supporting the variant">
##FORMAT=<ID=SC,Number=1,Type=Integer,Description="Number of soft-clipped alignments supporting the variant">
##FORMAT=<ID=BE,Number=1,Type=Integer,Description="Block edge metric">
##FORMAT=<ID=COV,Number=1,Type=Float,Description="Maximum read coverage +/- 10kb around break site at A or B">
##FORMAT=<ID=LNK,Number=1,Type=Integer,Description="Contig A and contig B overlap">
##FORMAT=<ID=NEIGH,Number=1,Type=Integer,Description="Number of other beak points within 100 bp or break sites">
##FORMAT=<ID=RB,Number=1,Type=Integer,Description="Number of reference bases in contigs">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Number of reads on plus strand">
##FORMAT=<ID=MS,Number=1,Type=Integer,Description="Number of reads on minus strand">{}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

# ##INFO=<ID=MPROB,Number=1,Type=Float,Description="Median probability of event across samples">
# ##FORMAT=<ID=PROB,Number=1,Type=Float,Description="Probability of event">

    outfile.write(HEADER.format(contig_names) + "\t" + "\t".join(names) + "\n")

    if show_names:
        click.echo("Input samples: {}".format(str(list(names))), err=True)

    version = pkg_resources.require("dysgu")[0].version

    if len(names) > 1:
        dm = df.sort_values(["partners"], ascending=False)
    else:
        dm = df

    seen_idx = set([])

    for col in ['raw_reads_10kb', 'DP', 'DN', 'DApri', 'DAsupp', 'NMpri', 'NMsupp', 'MAPQpri', 'MAPQsupp']:
        dm[col] = dm[col].round(2)

    for col in ['maxASsupp', 'neigh']:
        dm[col] = [int(i) for i in dm[col]]

    count = 0
    recs = []
    jobs = []

    add_kind = args["add_kind"] == "True"
    for idx, r in dm.iterrows():

        if idx in seen_idx:
            continue

        format_f, df_rows = gen_format_fields(r, df, names)

        if "partners" in r:
            seen_idx |= set(r["partners"])

        r_main = make_main_record(r, version, count, format_f, df_rows, add_kind)
        recs.append(r_main)
        count += 1

    for rec in sorted(recs, key=lambda x: (x[0], x[1])):
        outfile.write("\t".join(list(map(str, rec))) + "\n")