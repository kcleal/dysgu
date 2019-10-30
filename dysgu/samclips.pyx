#cython: language_level=3

"""
Utils to generate proper sam output and flag information
"""
from __future__ import absolute_import
import re
import click
from dysgu import io_funcs
from dysgu import samflags


def echo(*arg):
    click.echo(arg, err=True)


def set_tlen(out):

    pri_1 = out[0][1]
    pri_2 = out[1][1]

    p1_pos = int(pri_1[2])
    p2_pos = int(pri_2[2])

    flg1 = pri_1[0]
    flg2 = pri_2[0]

    assert flg1 & 64
    assert flg2 & 128

    # Use the end position of the alignment if read is on the reverse strand, or start pos if on the forward
    if flg1 & 16:
        aln_end1 = io_funcs.get_align_end_offset(pri_1[4])
        t1 = p1_pos + aln_end1
    else:
        t1 = p1_pos

    if flg2 & 16:
        aln_end2 = io_funcs.get_align_end_offset(pri_2[4])
        t2 = p2_pos + aln_end2
    else:
        t2 = p2_pos

    if t1 <= t2:
        tlen1 = t2 - t1  # Positive
        tlen2 = t1 - t2  # Negative
    else:
        tlen1 = t2 - t1  # Negative
        tlen2 = t1 - t2  # Positive

    pri_1[7] = str(tlen1)
    pri_2[7] = str(tlen2)

    out2 = [(out[0][0], pri_1, out[0][2]), (out[1][0], pri_2, out[1][2])]

    # Set tlen's of supplementary
    for sup_tuple in out[2:]:
        sup_tuple = list(sup_tuple)
        sup_flg = sup_tuple[1][0]
        sup_chrom = sup_tuple[1][1]
        sup_pos = int(sup_tuple[1][2])

        sup_end = io_funcs.get_align_end_offset(sup_tuple[1][4])
        if sup_flg & 16:  # If reverse strand, count to end
            sup_pos += sup_end

        if sup_flg & 64:  # First in pair, mate is second
            other_end = t2
            other_chrom = pri_2[1]
            other_flag = pri_2[0]
        else:
            other_end = t1
            other_chrom = pri_1[1]
            other_flag = pri_1[0]
        # This is a bit of a hack to make the TLEN identical to bwa
        # Make sure they are both on same chromosome
        if sup_chrom == other_chrom:
            if sup_pos < other_end:
                if bool(sup_flg & 16) != bool(other_flag & 16):  # Different strands
                    tlen = other_end - sup_pos
                else:
                    tlen = sup_pos - other_end
            else:
                if bool(sup_flg & 16) != bool(other_flag & 16):  # Different strands
                    tlen = other_end - sup_pos
                else:
                    tlen = sup_pos - other_end

            sup_tuple[1][7] = str(tlen)
        out2.append(tuple(sup_tuple))

    return out2


def set_mate_flag(a, b, max_d, read1_rev, read2_rev):

    if not a or not b:  # No alignment, mate unmapped?
        return False, False

    # Make sure chromosome of mate is properly set not "*"
    chrom_a, mate_a = a[2], a[5]
    chrom_b, mate_b = b[2], b[5]
    if chrom_a != mate_b:
        b[5] = chrom_a
    if chrom_b != mate_a:
        a[5] = chrom_b

    aflag = a[0]
    bflag = b[0]

    reverse_A = False
    reverse_B = False

    # If set as not primary, and has been aligned to reverse strand, and primary is mapped on forward
    # the sequence needs to be rev complement
    if aflag & 256:
        if (aflag & 16) and (not read1_rev):
            reverse_A = True
        elif (not aflag & 16) and read1_rev:
            reverse_A = True

    if bflag & 256:
        if (bflag & 16) and (not read2_rev):
            reverse_B = True
        elif (not bflag & 16) and read2_rev:
            reverse_B = True

    # Turn off proper pair flag, might be erroneously set
    aflag = samflags.set_bit(aflag, 1, 0)
    bflag = samflags.set_bit(bflag, 1, 0)

    # Turn off supplementary pair flag
    aflag = samflags.set_bit(aflag, 11, 0)
    bflag = samflags.set_bit(bflag, 11, 0)

    # Set paired
    aflag = samflags.set_bit(aflag, 0, 1)
    bflag = samflags.set_bit(bflag, 0, 1)

    # Set first and second in pair, in case not set
    aflag = samflags.set_bit(aflag, 6, 1)
    bflag = samflags.set_bit(bflag, 7, 1)

    # Turn off any mate reverse flags, these should be reset
    aflag = samflags.set_bit(aflag, 5, 0)
    bflag = samflags.set_bit(bflag, 5, 0)

    # If either read is unmapped
    if aflag & 4:
        bflag = samflags.set_bit(bflag, 3, 1)  # Position 3, change to 1
    if bflag & 4:
        aflag = samflags.set_bit(aflag, 3, 1)

    # If either read on reverse strand
    if aflag & 16:
        bflag = samflags.set_bit(bflag, 5, 1)
    if bflag & 16:
        aflag = samflags.set_bit(aflag, 5, 1)

    # Set unmapped
    arname = a[1]
    apos = a[2]
    if apos == "0":  # -1 means unmapped
        aflag = samflags.set_bit(aflag, 2, 1)
        bflag = samflags.set_bit(bflag, 8, 1)

    brname = b[1]
    bpos = b[2]
    if b[2] == "0":
        bflag = samflags.set_bit(bflag, 2, 1)
        aflag = samflags.set_bit(aflag, 8, 1)

    # Set RNEXT and PNEXT
    a[5] = brname
    a[6] = bpos

    b[5] = arname
    b[6] = apos

    if not (apos == "-1" or bpos == "-1"):

        if arname == brname:
            # Set TLEN
            p1, p2 = int(apos), int(bpos)

            # Set proper-pair flag
            if (aflag & 16 and not bflag & 16) or (not aflag & 16 and bflag & 16):  # Not on same strand

                if abs(p1 - p2) < max_d:
                    # Check for FR or RF orientation
                    if (p1 < p2 and (not aflag & 16) and (bflag & 16)) or (p2 <= p1 and (not bflag & 16) and (aflag & 16)):
                        aflag = samflags.set_bit(aflag, 1, 1)
                        bflag = samflags.set_bit(bflag, 1, 1)

                        # If proper pair, sometimes the mate-reverse-strand flag is set
                        # this subsequently means the sequence should be reverse complemented!
                        if aflag & 16 and not bflag & 32:
                            # Mate-reverse strand not set
                            bflag = samflags.set_bit(bflag, 5, 1)
                            # reverse_B = True

                        if not aflag & 16 and bflag & 32:
                            # Mate-reverse should'nt be set
                            bflag = samflags.set_bit(bflag, 5, 0)
                            reverse_A = True

                        if bflag & 16 and not aflag & 32:
                            # Mate-reverse strand not set
                            aflag = samflags.set_bit(aflag, 5, 1)
                            # reverse_A = True

                        if not bflag & 16 and aflag & 32:
                            # Mate-revsere should'nt be set
                            aflag = samflags.set_bit(aflag, 5, 0)
                            reverse_B = True

    a[0] = aflag
    b[0] = bflag
    return reverse_A, reverse_B


def set_supp_flags(sup, pri, ori_primary_reversed, primary_will_be_reversed):

    # Set paired
    supflag = sup[0]
    priflag = pri[0]

    # Set paired and supplementary flag
    if not supflag & 1:
        supflag = samflags.set_bit(supflag, 0, 1)
    if not supflag & 2048:
        supflag = samflags.set_bit(supflag, 11, 1)

    # If primary is on reverse strand, set the mate reverse strand tag
    if priflag & 16 and not supflag & 32:
        supflag = samflags.set_bit(supflag, 5, 1)
    # If primary is on forward srand, turn off mate rev strand
    if not priflag & 16 and supflag & 32:
        supflag = samflags.set_bit(supflag, 5, 0)

    # Turn off not-primary-alignment
    if supflag & 256:
        supflag = samflags.set_bit(supflag, 8, 0)

    rev_sup = False
    if ori_primary_reversed:
        if not supflag & 16:  # Read on forward strand
            rev_sup = True

    elif supflag & 16:  # Read on reverse strand
        if not ori_primary_reversed:
            rev_sup = True

    elif not supflag & 16:  # Read on forward strand
        if primary_will_be_reversed and not priflag & 16:  # Primary will end up on forward
            rev_sup = True  # Old primary on reverse, so needs rev comp

    sup[0] = supflag
    sup[5] = pri[1]
    sup[6] = pri[2]

    return rev_sup


def add_sequence_back(item, reverse_me, template):
    # item is the alignment
    flag = item[0]
    c = re.split(r'(\d+)', item[4])[1:]  # Drop leading empty string
    start = 0

    if flag & 64:  # Read1
        seq = template["read1_seq"]
        q = template["read1_q"]
    elif flag & 128:
        seq = template["read2_seq"]
        q = template["read2_q"]
    else:
        seq = template["read1_seq"]  # Unpaired
        q = template["read1_q"]

    cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "DH"])

    if len(seq) != cigar_length:
        if template["replace_hard"] and q != "*":
            # Sometimes current read had a hard-clip in cigar, but the primary read was not soft clipped
            cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "D"])
            echo(c, cigar_length, len(seq), len(template["read1_q"]), template['name'], template["read1_q"])
            if len(seq) != cigar_length:
                return item  # Cigar length is not set properly by mapper
            # If this is true, reset the Hard-clips with Soft-clips
            item[4] = item[4].replace("H", "S")
        return item

    # Occasionally the H is missing, means its impossible to add sequence back in

    total_cigar_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1]])
    if (flag & 64 and len(template["read1_seq"]) > total_cigar_length) or \
            (flag & 128 and len(template["read2_seq"]) > total_cigar_length):
        return item

    if flag & 64 and template["read1_seq"]:
        name = "read1"
        if template["fq_read1_seq"] != 0:
            end = len(template["fq_read1_seq"])
        else:
            end = len(template["read1_seq"])

    elif flag & 128 and template["read2_seq"]:
        name = "read2"
        if template["fq_read2_seq"] != 0:
            end = len(template["fq_read2_seq"])
        else:
            end = len(template["read2_seq"])
    else:
        return item  # read sequence is None or bad flag

    # Try and replace H with S
    if c[1] == "H" or c[-1] == "H":
        # Replace hard with soft-clips
        non_hardclipped_length = sum([int(c[i]) for i in range(0, len(c), 2) if c[i + 1] not in "D"])
        if non_hardclipped_length == end and template["replace_hard"]:
            item[4] = item[4].replace("H", "S")
            cigar_length = non_hardclipped_length
        else:
            # Remove seq
            if c[1] == "H":
                start += int(c[0])
            if c[-1] == "H":
                end -= int(c[-2])

    # Might need to collect from the reverse direction; swap end and start
    if flag & 256 or flag & 2048:
        if flag & 64 and template["read1_reverse"] != bool(flag & 16):
            # Different strand to primary, count from end
            new_end = template["read1_length"] - start
            new_start = template["read1_length"] - end
            start = new_start
            end = new_end

        elif flag & 128 and (template["read2_reverse"] != bool(flag & 16)):
            new_end = template["read2_length"] - start
            new_start = template["read2_length"] - end
            start = new_start
            end = new_end

    # Try and use the primary sequence to replace hard-clips
    if item[9] == "*" or len(item[9]) < abs(end - start) or len(item[9]) == 0:
        if template["replace_hard"] and template["fq_%s_q" % name]:
            key = "fq_"
        else:
            key = ""
        s = template["%s%s_seq" % (key, name)][start:end]
        q = template["%s%s_q" % (key, name)][start:end]

        if len(s) == cigar_length:
            if reverse_me:
                item[8] = io_funcs.reverse_complement(s, len(s))
                item[9] = q[::-1]
            else:
                item[8] = s
                item[9] = q

    # Try and use the supplied fq file to replace the sequence
    elif template["fq_%s_q" % name] != 0 and len(template["fq_%s_q" % name]) > len(item[9]):
        if item[9] in template["fq_%s_q" % name]:
            item[8] = template["fq_%s_seq" % name][start:end]
            item[9] = template["fq_%s_q" % name][start:end]

        elif item[9] in template["fq_%s_q" % name][::-1]:
            sqn = "fq_%s_seq" % name
            s = io_funcs.reverse_complement(template[sqn], len(template[sqn]))[start:end]
            q = template["fq_%s_q" % name][::-1][start:end]
            if len(s) == cigar_length:
                item[8] = s
                item[9] = q

        else:
            echo("---")
            echo(item[9], flag)
            echo(name)
            echo(template["read1_q"])
            echo(template["read2_q"])
            echo(item)
            raise ValueError

    if len(item[8]) != cigar_length:
        echo(len(item[8]), cigar_length, len(item[9]), start, end)
        echo(template)

    assert len(item[8]) == cigar_length
    return item


def replace_sa_tags(alns):

    if any([i[0] == "sup" for i in alns]):
        sa_tags = {}  # Read1: tag, might be multiple split alignments
        alns2 = []
        for i, j, k in alns:
            # Remove any SA tags in alignment, might be wrong
            j = [item for idx, item in enumerate(j) if idx <= 9 or (idx > 9 and item[:2] != "SA")]
            flag = j[0]
            mapq = j[3]
            nm = 0
            chrom = j[1]
            pos = j[2]
            for tg in j[10:]:
                if tg[:2] == "NM":
                    nm = tg[5:]
                    break

            strand = "-" if flag & 16 else "+"
            cigar = j[4]
            sa = "%s,%s,%s,%s,%s,%s,%s" % (chrom, pos, strand, cigar, j[0], mapq, nm)
            key = (flag & 64, 1 if flag & 2048 else 0)
            if key in sa_tags:
                sa_tags[key] += ";" + sa
            else:
                sa_tags[key] = sa
            alns2.append([i, j, k])

        # Now add back in
        out = []
        for i, j, k in alns2:
            flag = j[0]
            key = (flag & 64, 0 if flag & 2048 else 1)
            if key in sa_tags:
                j.insert(14, "SA:Z:" + sa_tags[key])
            out.append((i, j, k))
        return out
    else:
        # Might need to remove SA tags
        return [(i, [item for idx, item in enumerate(j) if idx <= 9 or (idx > 9 and item[:2] != "SA")], ii) for i, j, ii in alns]


def fixsam(template):

    sam = [template['inputdata'][i] for i in template['rows']]
    max_d = template['max_d']

    # Todo make sure all read-pairs have a mapping, otherwise write an unmapped

    paired = False if template["read2_length"] is 0 else True
    score_mat = template['score_mat']

    out = []
    primary1 = 0
    primary2 = 0
    rev_A = False
    rev_B = False
    for l in sam:

        l[0] = int(l[0])  # Convert flag to int
        t = score_mat
        strand = "-1" if l[0] & 16 else "1"
        rid = str(2 if l[0] & 128 else 1)
        key = "{}-{}-{}-{}".format(l[1], l[2], strand, rid)

        if len(t[key]) > 2:
            # Prevent bug where two identical alignments possible
            aln_info_0, aln_info_1 = t[key].pop(0), t[key].pop(0)  # Remove first two items from list
        else:
            aln_info_0, aln_info_1 = t[key]

        xs = int(aln_info_1)
        if l[0] & 2048:
            os = "DS:i:1"  # refers to "originally supplementary"
        else:
            os = "DS:i:0"
        l += [
              "DA:i:" + str(xs),
              "DP:Z:" + str(round(t["dis_to_next_path"], 0)),
              "DN:Z:" + str(round(t["dis_to_normal"], 2)),
              "PS:Z:" + str(round(t["path_score"], 2)),
              "NP:Z:" + str(round(t["normal_pairings"], 1)),
              os
              ]

        if aln_info_0:
            if rid == "1":
                primary1 = l
            else:
                primary2 = l
        else:
            out.append(['sup', l, False])  # Supplementary, False to decide if rev comp

    if primary1 is 0 or primary2 is 0 and template["paired_end"]:
        return []  # Todo deal with unmapped read or unpaired

    if paired and template["paired_end"]:
        rev_A, rev_B = set_mate_flag(primary1, primary2, max_d, template["read1_reverse"], template["read2_reverse"])

        # Check if supplementary needs reverse complementing
        for i in range(len(out)):
            if out[i][1][0] & 64:  # First in pair  Note primary2 and primary1 were switched
                revsup = set_supp_flags(out[i][1], primary2, template["read1_reverse"], rev_A)
            else:
                revsup = set_supp_flags(out[i][1], primary1, template["read2_reverse"], rev_B)

            if revsup:
                out[i][2] = True

    if template["paired_end"]:
        out = [('pri', primary1, rev_A), ('pri', primary2, rev_B)] + out
        out = set_tlen(out)

    else:
        out = [('pri', primary1, rev_A)] + out

    # Add read seq info back in if necessary, before reverse complementing. Check for hard clips and clip as necessary

    for a_type, aln, reverse_me in out:

        if aln:  # None here means no alignment for primary2

            # Do for all supplementary
            if len(aln[8]) <= 1 or "H" in aln[4] or aln[0] & 2048:  # Sequence is set as "*", needs adding back in

                aln = add_sequence_back(aln, reverse_me, template)

            elif reverse_me:
                aln[8] = io_funcs.reverse_complement(str(aln[8]), len(aln[8]))
                aln[9] = aln[9][::-1]

            # Turn off not primary here
            aln[0] = samflags.set_bit(aln[0], 8, 0)

    out = replace_sa_tags(out)

    # Convert flag back to string
    for j in range(len(out)):
        out[j][1][0] = str(out[j][1][0])

    return [i[1] for i in out if i[1] != 0]

