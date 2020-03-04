# cython: language_level=3

from __future__ import absolute_import
from collections import Counter, defaultdict
import click
import numpy as np
from dysgu import data_io, assembler, graph, coverage
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


def echo(*args):
    click.echo(args, err=True)


cdef dict count_attributes(list reads1, list reads2, int min_support):

    if len(reads1) == 0:
        raise ValueError("No reads in set")

    r = {"pe": 0, "supp": 0, "sc": 0, "DP": [], "DApri": [], "DN": [], "NMpri": [], "NP": 0, "DAsupp": [], "NMsupp": [],
         "maxASsupp": [], "MAPQpri": [], "MAPQsupp": [], "plus": 0, "minus": 0}
    paired_end = set([])
    seen = set([])

    cdef int flag
    cdef list reads
    for reads in (reads1, reads2):
        for a in reads:

            qname = a.qname
            if qname not in seen:  # Add these once for each pair, its common across all alignments of template
                if a.has_tag("DP"):
                    r["DP"].append(float(a.get_tag("DP")))
                if a.has_tag("DN"):
                    r["DN"].append(float(a.get_tag("DN")))
                if a.has_tag("NP"):
                    if float(a.get_tag("NP")) == 1:
                        r["NP"] += 1
                seen.add(qname)

            flag = a.flag
            if flag & 2048:  # Supplementary
                r["supp"] += 1
                r["MAPQsupp"].append(a.mapq)
                if a.has_tag("DA"):
                    r["DAsupp"].append(float(a.get_tag("DA")))
                if a.has_tag("NM"):
                    r["NMsupp"].append(float(a.get_tag("NM")))
                if a.has_tag("AS"):
                    r["maxASsupp"].append(float(a.get_tag("AS")))

            elif not flag & 256:  # Primary reads
                r["MAPQpri"].append(a.mapq)
                if qname in paired_end:  # If two primary reads from same pair
                    r["pe"] += 1
                else:
                    paired_end.add(qname)
                if a.has_tag("DA"):
                    r["DApri"].append(float(a.get_tag("DA")))
                if a.has_tag("NM"):
                    r["NMpri"].append(float(a.get_tag("NM")))

            if flag & 16:
                r["minus"] += 1
            else:
                r["plus"] += 1

            ct = a.cigartuples
            if ct[0][0] == 4 or ct[-1][0] == 4:
                r["sc"] += 1

    if r["pe"] + r["supp"] < min_support:
        return {}

    cdef str k
    for k in ("DP", "DApri", "DN", "NMpri", "DAsupp", "NMsupp", "MAPQpri", "MAPQsupp"):
        if len(r[k]) > 0:
            r[k] = np.mean(r[k])
        else:
            r[k] = 0

    if len(r["maxASsupp"]) > 0:
        r["maxASsupp"] = int(max(r["maxASsupp"]))
    else:
        r["maxASsupp"] = 0

    return r


cdef dict fetch_reads(dict data, d, bam):

    input_reads = data["reads"]
    dta_n2n = data["n2n"]
    if len(dta_n2n) > 0:  # Might need to collect reads from file
        n2n = {}  # Subset
        for block in d.values():
            for v in block:
                if v in dta_n2n:
                    n2n[v] = dta_n2n[v]  # Needs collecting

        input_reads.update(graph.get_reads(bam, n2n))

    return input_reads


def guess_informative_pair(aligns):

    if len(aligns) == 2:
        # Make sure aligns map different break points
        a = aligns[0]
        b = aligns[1]
        if abs(a.pos - b.pos) < 25 or abs(a.reference_end - b.reference_end) < 25:
            return
        if a.pos < b.pos:  # a and b will be same on same chrom
            return a, b
        return b, a

    pri_first = None
    sup_first = None
    pri_second = None
    sup_second = None
    for i in aligns:
        if not i.flag & 2304:  # Not pri, supplementary --> is primary
            if i.flag & 64:
                pri_first = i
            else:
                pri_second = i
        elif i.flag & 2048:  # Supplementary
            if i.flag & 64:
                sup_first = i
            else:
                sup_second = i
    a = None
    b = None
    if pri_first and sup_first:
        a = pri_first
        b = sup_first
    elif pri_second and sup_second:
        a = pri_second
        b = sup_second
    elif pri_first and pri_second:
        a = pri_first
        b = pri_second
    if a is None:
        return None
    if a.pos < b.pos:
        return a, b
    return b, a


cdef dict single(infile, dict data, int insert_size, int insert_stdev, int clip_length, int min_support,
                 int to_assemble=1):
    # Make sure at least one read is worth calling
    cdef int min_distance = insert_size + (2*insert_stdev)
    cdef int n_templates = len(set([i.qname for i in data["reads"].values()]))
    # if n_templates < min_support:
    #     return {}

    if n_templates == 1:
        if not any((not i.flag & 2) or (i.rname != i.rnext) or
                   (abs(i.tlen) > min_distance) for i in data["reads"].values()):
            return {}

    # Single's can sometimes be seperated into two groups --> send one_edge
    # Otherwise, infer the other breakpoint from a single group

    # Groupby template name
    precise_a = []
    precise_b = []
    informative = []
    u_reads = []
    v_reads = []
    tmp = defaultdict(list)
    for align in data["reads"].values():
        tmp[align.qname].append(align)

    cdef AlignmentItem v_item, itm
    cdef int left_clip_a, right_clip_a, left_clip_b, right_clip_b

    for temp_name, alignments in tmp.items():
        l_align = list(alignments)
        if len(l_align) > 1:
            pair = guess_informative_pair(l_align)
            if pair:
                u_reads.append(pair[0])
                v_reads.append(pair[1])

                a, b = pair
                a_ct = a.cigartuples
                b_ct = b.cigartuples
                # Soft-clips for the chosen pair, plus template start of alignment
                left_clip_a, right_clip_a, left_clip_b, right_clip_b = mask_soft_clips(a.flag, b.flag, a_ct, b_ct)

                v_item = AlignmentItem(a.rname,
                                       b.rname,
                                       int(not a.flag & 3328),
                                       int(not b.flag & 3328),
                                       1 if a.flag & 64 else 2,
                                       1 if b.flag & 64 else 2,
                                       a.pos,
                                       a.reference_end,
                                       b.pos,
                                       b.reference_end,
                                       3 if a.flag & 16 == 0 else 5,
                                       3 if b.flag & 16 == 0 else 5,
                                       left_clip_a,
                                       right_clip_a,
                                       left_clip_b,
                                       right_clip_b,
                                       )
                
                if v_item.left_clipA and v_item.right_clipA:
                    if a_ct[0][1] >= a_ct[-1][1]:
                        v_item.right_clipA = 0
                    else:
                        v_item.left_clipA = 0

                if v_item.left_clipB and v_item.right_clipB:
                    if b_ct[0][1] >= b_ct[-1][1]:
                        v_item.right_clipB = 0
                    else:
                        v_item.left_clipB = 0

                classify_d(v_item)

                if v_item.breakA_precise:
                    precise_a.append(v_item.breakA)
                if v_item.breakB_precise:
                    precise_b.append(v_item.breakB)

                informative.append(v_item)

    if not informative:
        return {}

    attrs = count_attributes(u_reads, v_reads, min_support)

    if not attrs:
        return {}

    call_informative = Counter([(itm.svtype, itm.join_type) for itm in informative]).most_common()

    cdef str svtype, jointype

    svtype, jointype = call_informative[0][0]
    info = make_call(informative, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev)

    if not info:
        return {}

    info.update(attrs)

    as1 = None
    as2 = None

    if info["preciseA"]:
        as1 = assembler.base_assemble(u_reads)

    if info["preciseB"]:
        as2 = assembler.base_assemble(v_reads)

    info["linked"] = 0
    info["block_edge"] = 0
    info["contig"] = None
    info["contig2"] = None
    rbases = 0
    if as1 is not None and "contig" in as1:
        info["contig"] = as1["contig"]
        rbases += as1["ref_bases"]

    if as2 is not None and "contig" in as2:
        info["contig2"] = as2["contig"]
        rbases += as2["ref_bases"]
    info["ref_bases"] = rbases

    return info


def informative_pair(u, v):

    pri_u = None
    pri_v = None
    for i in u:
        ri_flag = i.flag & 64
        if not i.flag & 2304:  # Not pri, supplementary --> is primary
            pri_u = i
        for j in v:
            if j.flag & 64 == ri_flag:  # Same read
                # Same read, primary + supp, or supp + supp
                return i, j

            if not j.flag & 2304:  # Is primary
                pri_v = j

            # else Different read and sup, e.g. read 1 primary + read 2 supplementary (not informative)

    if pri_u is not None and pri_v is not None:
        return pri_u, pri_v


cdef void two_primary(AlignmentItem v):

    if v.posA < v.posB or (v.posA == v.posB and v.endA < v.endB):

        if v.strandA == 3 and v.strandB == 5:  # DEL type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "DEL"
            v.join_type = "3to5"

        elif v.strandA == 5 and v.strandB == 3:  # DUP type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

        elif v.strandA == v.strandB:  # INV type
            # Break to left
            if v.strandA == 5:

                if not (v.left_clipA or v.left_clipB) and (v.right_clipA or v.right_clipB):
                    v.breakA = v.endA
                    if v.right_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.endB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "3to3"

                else:
                    v.breakA = v.posA
                    if v.left_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.posB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "5to5"

            elif v.strandA == 3:

                if not (v.right_clipA or v.right_clipB) and (v.left_clipA or v.left_clipB):
                    v.breakA = v.posA
                    if v.right_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.posB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "5to5"

                else:

                    v.breakA = v.endA
                    if v.left_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.endB
                    if v.right_clipB:
                        v.breakB_precise = 1

                    v.svtype = "INV"
                    v.join_type = "3to3"

            # Break INV:DUP
            else:
            # elif v.left_clipA and v.right_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV:DUP"
                v.join_type = "3to5"

    else:  # B < A

        if v.strandA == 5 and v.strandB == 3:  # DEL type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "DEL"
            v.join_type = "3to5"

        elif v.strandB == 5 and v.strandA == 3:  # DUP type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

            # INV type
        elif v.strandA == 5:

            if not (v.left_clipA or v.left_clipB) and (v.right_clipA or v.right_clipB):
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "3to3"

            else:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.posB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "5to5"

        elif v.strandA == 3:

            if not (v.right_clipA or v.right_clipB) and (v.left_clipA or v.left_clipB):
                v.breakA = v.posA
                if v.right_clipA:
                    v.breakA_precise = 1

                v.breakB = v.posB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "5to5"

            else:

                v.breakA = v.endA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "3to3"


cdef void same_read(AlignmentItem v):
    # A is first
    if v.posA < v.posB or (v.posA == v.posB and v.endA <= v.endB):

        if v.strandA == v.strandB:  # Same for 3 and 5 strand reads

            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested DUP
                if v.left_clipA:
                    v.breakA = v.posA
                    v.breakA_precise = 1
                elif v.right_clipA:
                    v.breakA = v.endA
                    v.breakA_precise = 1
                else:
                    v.breakA = v.endA

                if v.left_clipB:
                    v.breakB = v.posB
                    v.breakB_precise = 1
                elif v.right_clipB:
                    v.breakB = v.endB
                    v.breakB_precise = 1
                else:
                    v.breakB = v.endB
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif not v.left_clipA and not v.right_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "DEL"
                v.join_type = "3to5"

            elif not v.right_clipA and not v.left_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif not v.left_clipA and not v.left_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            else:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"

        elif v.strandA == 5 and v.strandB == 3:

            # Break right
            if not v.left_clipA and not v.left_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "3to3"

            # Break left
            elif not v.right_clipA and not v.right_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1

                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1

                v.svtype = "INV"
                v.join_type = "5to5"

            elif is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Inverted duplication
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV:DUP"
                v.join_type = "5to3"

            elif v.right_clipA and v.left_clipB:
                v.breakA = v.endA
                v.breakB = v.posB
                v.svtype = "INV"
                v.join_type = "3to5"

            else:
            # elif v.left_clipA and v.right_clipB:
                v.breakA = v.posA
                v.breakB = v.endB
                v.svtype = "INV"
                v.join_type = "5to3"

        else:  # INV type
            # Break left
            if v.left_clipA and v.left_clipB:
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                v.join_type = "5to5"
            # Break right
            elif v.right_clipA and v.right_clipB:
                v.breakA = v.endA
                v.breakB = v.endB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            else:  # Guess using pos only
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                if v.strandA == 5:
                    v.join_type = "5to5"
                else:
                    v.join_type = "3to3"

    # B is first
    else:
        if v.strandA == v.strandB:

            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested DUP
                if v.left_clipA:
                    v.breakA = v.posA
                    v.breakA_precise = 1
                elif v.right_clipA:
                    v.breakA = v.endA
                    v.breakA_precise = 1
                else:
                    v.breakA = v.endA

                if v.left_clipB:
                    v.breakB = v.posB
                    v.breakB_precise = 1
                elif v.right_clipB:
                    v.breakB = v.endB
                    v.breakB_precise = 1
                else:
                    v.breakB = v.endB
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif not v.left_clipB and not v.right_clipA:
                v.breakA = v.posA
                v.breakB = v.endB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DEL"
                v.join_type = "3to5"

            elif not v.right_clipB and not v.left_clipA:
                v.breakA = v.endA
                v.breakB = v.posB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            else:
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "BND"
                v.join_type = f"{v.strandA}to{v.strandB}"

        else:  # INV type
            # Break left
            if v.left_clipA and v.left_clipB:
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                v.join_type = "5to5"

            # Break right
            elif v.right_clipA and v.right_clipB:
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.breakA = v.endA
                v.breakB = v.endB
                v.svtype = "INV"
                v.join_type = "3to3"

            else:  # Guess using pos only
                v.breakA = v.posA
                v.breakB = v.posB
                v.svtype = "INV"
                if v.strandA == 5:
                    v.join_type = "5to5"
                else:
                    v.join_type = "3to3"


cdef void different_read(AlignmentItem v):

    if v.posA < v.posB or (v.posA == v.posB and v.endA <= v.endB):  # A is first

        if v.strandA == 3 and v.strandB == 5:  # DEL type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "DEL"
            v.join_type = "3to5"

        elif v.strandA == 5 and v.strandB == 3:  # DUP type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

        elif v.strandA == v.strandB:
            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested

                if v.strandA == 3:  # Both forward strand
                    v.breakA = v.endA
                    if v.right_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.endB
                    if v.right_clipB:
                        v.breakB_precise = 1
                    v.svtype = "INV"
                    v.join_type = "3to3"

                else:  # Both forward strand
                    v.breakA = v.posA
                    if v.left_clipA:
                        v.breakA_precise = 1

                    v.breakB = v.posB
                    if v.left_clipB:
                        v.breakB_precise = 1
                    v.svtype = "INV"
                    v.join_type = "5to5"

            elif not v.right_clipA and not v.right_clipB:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"

            elif not v.left_clipA and not v.left_clipB:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            elif v.right_clipA and v.left_clipB:
                v.breakA = v.endA
                v.breakA_precise = 1
                v.breakB = v.posB
                v.breakB_precise = 1
                v.svtype = "INV:DUP"
                v.join_type = "5to3"

            else:
            # elif v.left_clipA and v.right_clipB:
                v.breakA = v.posA
                v.breakA_precise = 1
                v.breakB = v.endB
                v.breakB_precise = 1
                v.svtype = "INV:DUP"
                v.join_type = "3to5"

    else:  # B is first; B <= A

        if v.strandA == 5 and v.strandB == 3:  # DEL type
            v.breakA = v.posA
            if v.left_clipA:
                v.breakA_precise = 1

            v.breakB = v.endB
            if v.right_clipB:
                v.breakB_precise = 1

            v.svtype = "DEL"
            v.join_type = "3to5"

        elif v.strandA == 3 and v.strandB == 5:  # DUP type
            v.breakA = v.endA
            if v.right_clipA:
                v.breakA_precise = 1

            v.breakB = v.posB
            if v.left_clipB:
                v.breakB_precise = 1

            v.svtype = "DUP"
            v.join_type = "5to3"

        elif v.strandA == v.strandB:  # INV type

            if is_overlapping(v.posA, v.endA, v.posB, v.endB):  # Nested DUP
                if v.left_clipA:
                    v.breakA = v.posA
                    v.breakA_precise = 1
                elif v.right_clipA:
                    v.breakA = v.endA
                    v.breakA_precise = 1
                else:
                    v.breakA = v.endA

                if v.left_clipB:
                    v.breakB = v.posB
                    v.breakB_precise = 1
                elif v.right_clipB:
                    v.breakB = v.endB
                    v.breakB_precise = 1
                else:
                    v.breakB = v.endB
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif v.right_clipA and v.left_clipB:
                v.breakA = v.endA
                v.breakB = v.posB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DEL"
                v.join_type = "3to5"

            elif v.left_clipA and v.right_clipB:
                v.breakA = v.posA
                v.breakB = v.endB
                v.breakA_precise = 1
                v.breakB_precise = 1
                v.svtype = "DUP"
                v.join_type = "5to3"

            elif not v.left_clipB and not v.left_clipA:

                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1

                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            elif not v.right_clipB and not v.right_clipA:

                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1

                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"

            elif v.strandA == 3:
                v.breakA = v.endA
                if v.right_clipA:
                    v.breakA_precise = 1
                v.breakB = v.endB
                if v.right_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "3to3"

            else:
                v.breakA = v.posA
                if v.left_clipA:
                    v.breakA_precise = 1
                v.breakB = v.posB
                if v.left_clipB:
                    v.breakB_precise = 1
                v.svtype = "INV"
                v.join_type = "5to5"


cdef void translocation(AlignmentItem v):

    v.svtype = "TRA"

    if v.left_clipA:
        v.breakA = v.posA
        v.breakA_precise = 1
    elif v.right_clipA:
        v.breakA = v.endA
        v.breakA_precise = 1
    else:
        v.breakA = v.posA
        v.breakA_precise = 1

    if v.left_clipB:
        v.breakB = v.posB
        v.breakB_precise = 1
    elif v.right_clipB:
        v.breakB = v.endB
        v.breakB_precise = 1
    else:
        v.breakB = v.posB
        v.breakB_precise = 1
    v.join_type = f"{v.strandA}to{v.strandB}"


cdef int is_overlapping(int x1, int x2, int y1, int y2) nogil:
    return int(max(x1, y1) <= min(x2, y2))


cdef void classify_d(AlignmentItem v):

    v.breakA_precise = 0
    v.breakB_precise = 0

    if v.chrA != v.chrB:
        translocation(v)

    else:  # Intra-chromosomal
        # Find join type first, different for pri-sup, pri-pri
        if v.priA and v.priB:  # Both primary
            two_primary(v)

        else:  # One is a supplementary
            if v.rA == v.rB:
                same_read(v)
            else:
                different_read(v)


cdef tuple break_ops(list positions, list precise, int limit, float median_pos):

    ops = []
    cdef int i
    for i in positions:
        ops.append((i, 1, 0))
        ops.append((i + limit, -1, 0))
    cdef int v = 1
    if limit < 0:
        v = -1
    for i in precise:
        ops.append((i, 1, 1))
        if v == -1:
            ops.append((i + v + 1, -1, 1))
        else:
            ops.append((i + v, -1, 1))

    ops = sorted(ops, reverse=limit < 0)
    cdef int cum_sum = 0
    cdef int max_sum = 0
    cdef int max_sum_i = 0
    cdef int idx
    for idx, item in enumerate(ops):
        cum_sum += item[1]
        if cum_sum > max_sum:
            max_sum = cum_sum
            max_sum_i = idx

    cdef int break_point, cipos95, is_precise
    break_point, _, is_precise = ops[max_sum_i]
    cipos95 = 0
    if not is_precise and not any(i == break_point for i in precise):
        # Calculate confidence interval around break
        cipos95 = int(abs(int(np.percentile(positions, [97.5])) - median_pos))

    return break_point, cipos95


cdef dict make_call(list informative, list breakA_precise, list breakB_precise, str svtype, str jointype,
                    int insert_size, int insert_stdev):
    # Inspired by mosdepth algorithm +1 for start -1 for end using intervals where break site could occur
    # Use insert size to set a limit on where break site could occur
    cdef int limit = insert_size + insert_stdev

    # get bulk call
    positionsA = [i.breakA for i in informative]
    positionsB = [i.breakB for i in informative]
    cdef float median_A = np.median(positionsA)
    cdef float median_B = np.median(positionsB)

    cdef int main_A_break = 0
    cdef int main_B_break = 0
    cdef int cipos95A, cipos95B

    if svtype == "DEL":
        if median_A < median_B:
            main_A_break, cipos95A = break_ops(positionsA, breakA_precise, limit, median_A)
            main_B_break, cipos95B = break_ops(positionsB, breakB_precise, -limit, median_B)
        else:
            main_B_break, cipos95B = break_ops(positionsB, breakB_precise, limit, median_B)
            main_A_break, cipos95A = break_ops(positionsA, breakA_precise, -limit, median_A)

    elif svtype == "DUP":
        if median_A < median_B:
            main_A_break, cipos95A = break_ops(positionsA, breakA_precise, -limit, median_A)
            main_B_break, cipos95B = break_ops(positionsB, breakB_precise, limit, median_B)
        else:
            main_A_break, cipos95A = break_ops(positionsA, breakA_precise, limit, median_A)
            main_B_break, cipos95B = break_ops(positionsB, breakB_precise, -limit, median_B)

    # INV and TRA types
    elif jointype == "3to3":
        main_A_break, cipos95A = break_ops(positionsA, breakA_precise, limit, median_A)
        main_B_break, cipos95B = break_ops(positionsB, breakB_precise, limit, median_B)
    elif jointype == "5to5":
        main_A_break, cipos95A = break_ops(positionsA, breakA_precise, -limit, median_A)
        main_B_break, cipos95B = break_ops(positionsB, breakB_precise, -limit, median_B)
    # Non-canonical
    elif jointype == "3to5":
        main_A_break, cipos95A = break_ops(positionsA, breakA_precise, limit, median_A)
        main_B_break, cipos95B = break_ops(positionsB, breakB_precise, -limit, median_B)
    else:
        main_A_break, cipos95A = break_ops(positionsA, breakA_precise, -limit, median_A)
        main_B_break, cipos95B = break_ops(positionsB, breakB_precise, limit, median_B)

    return {"svtype": svtype, "join_type": jointype, "chrA": informative[0].chrA, "chrB": informative[0].chrB,
            "cipos95A": cipos95A, "cipos95B": cipos95B, "posA": main_A_break, "posB": main_B_break,
            "preciseA": True if cipos95A == 0 else False, "preciseB": True if cipos95B == 0 else False}


cdef tuple mask_soft_clips(int aflag, int bflag, list a_ct, list b_ct):

    # Find out which soft clip pairs are compatible with the chosen read pair

    cdef int left_clipA = 1 if (a_ct[0][0] == 4 or a_ct[0][0] == 5) else 0
    cdef int right_clipA = 1 if (a_ct[-1][0] == 4 or a_ct[-1][0] == 5) else 0
    cdef int left_clipB = 1 if (b_ct[0][0] == 4 or b_ct[0][0] == 5) else 0
    cdef int right_clipB = 1 if (b_ct[-1][0] == 4 or b_ct[-1][0] == 5) else 0

    cdef int a_template_start = 0
    cdef int b_template_start = 0

    if aflag & 64 == bflag & 64:  # Same read

        if (left_clipB and right_clipB) or (left_clipA and right_clipA):  # One read has more than one soft-clip

            if aflag & 16:  # A is reverse, convert to forward
                if a_ct[-1][0] != 0:
                    a_template_start = a_ct[-1][1]
            else:
                if a_ct[0][0] != 0:
                    a_template_start = a_ct[0][1]

            if bflag & 16:  # B is reverse
                if b_ct[-1][0] != 0:
                    b_template_start = b_ct[-1][1]
            else:
                if b_ct[0][0] != 0:
                    b_template_start = b_ct[0][1]

            if left_clipB and right_clipB:  # Choose one soft-clip for B

                if b_template_start < a_template_start:
                    if not bflag & 16:  # B on forward strand
                        left_clipB = 0
                    else:
                        right_clipB = 0

                else:
                    if not bflag & 16:  # B on reverse strand
                        right_clipB = 0
                    else:
                        left_clipB = 0

            if left_clipA and right_clipA:  # Choose one soft-clip for A

                if a_template_start < b_template_start:
                    if not aflag & 16:  # A on forward strand
                        left_clipA = 0
                    else:
                        right_clipA = 0

                else:
                    if not aflag & 16:  # A on reverse strand
                        right_clipA = 0
                    else:
                        left_clipA = 0

    else:  # Different reads choose longest if more than one soft-clip
        if left_clipB and right_clipB:
            if b_ct[0][1] > b_ct[-1][1]:
                right_clipB = 0
            else:
                left_clipB = 0

        if left_clipA and right_clipA:
            if a_ct[0][1] > a_ct[-1][1]:
                right_clipA = 0
            else:
                left_clipA = 0

    return left_clipA, right_clipA, left_clipB, right_clipB


cdef class AlignmentItem:
    cdef public int chrA, chrB, priA, priB, rA, rB, posA, endA, posB, endB, strandA, strandB, left_clipA, right_clipA,\
        left_clipB, right_clipB, breakA_precise, breakB_precise, breakA, breakB
    cdef public str svtype, join_type
    def __cinit__(self, int chrA, int chrB, int priA, int priB, int rA, int rB, int posA, int endA, int posB, int endB,
                  int strandA, int strandB, int left_clipA, int right_clipA, int left_clipB, int right_clipB,):
        self.chrA = chrA
        self.chrB = chrB
        self.priA = priA
        self.priB = priB
        self.rA = rA
        self.rB = rB
        self.posA = posA
        self.endA = endA
        self.posB = posB
        self.endB = endB
        self.strandA = strandA
        self.strandB = strandB
        self.left_clipA = left_clipA
        self.right_clipA = right_clipA
        self.left_clipB = left_clipB
        self.right_clipB = right_clipB

        #
        self.breakA_precise = 0
        self.breakB_precise = 0
        self.breakA = -1
        self.breakB = -1
        self.svtype = ""
        self.join_type = ""


cdef dict call_from_reads(u_reads, v_reads, int insert_size, int insert_stdev):

    grp_u = defaultdict(list)
    grp_v = defaultdict(list)

    for r in u_reads:
        grp_u[r.qname].append(r)
    for r in v_reads:
        grp_v[r.qname].append(r)

    informative = []
    precise_a = []
    precise_b = []

    cdef AlignmentItem v_item, i
    cdef int left_clip_a, right_clip_a, left_clip_b, right_clip_b

    for qname in [k for k in grp_u if k in grp_v]:  # Qname fount on both sides

        u = grp_u[qname]
        v = grp_v[qname]

        if len(u) == 1 and len(v) == 1:
            pair = (u[0], v[0])
        else:
            pair = informative_pair(u, v)

        if pair:

            a, b = pair
            a_ct = a.cigartuples
            b_ct = b.cigartuples
            # Soft-clips for the chosen pair, plus template start of alignment
            left_clip_a, right_clip_a, left_clip_b, right_clip_b = mask_soft_clips(a.flag, b.flag, a_ct, b_ct)

            v_item = AlignmentItem(a.rname,
                                   b.rname,
                                   int(not a.flag & 3328),
                                   int(not b.flag & 3328),
                                   1 if a.flag & 64 else 2,
                                   1 if b.flag & 64 else 2,
                                   a.pos,
                                   a.reference_end,
                                   b.pos,
                                   b.reference_end,
                                   3 if a.flag & 16 == 0 else 5,
                                   3 if b.flag & 16 == 0 else 5,
                                   left_clip_a,
                                   right_clip_a,
                                   left_clip_b,
                                   right_clip_b,
                                   )

            if v_item.left_clipA and v_item.right_clipA:
                if a_ct[0][1] >= a_ct[-1][1]:
                    v_item.right_clipA = 0
                else:
                    v_item.left_clipA = 0

            if v_item.left_clipB and v_item.right_clipB:
                if b_ct[0][1] >= b_ct[-1][1]:
                    v_item.right_clipB = 0
                else:
                    v_item.left_clipB = 0

            classify_d(v_item)

            if v_item.breakA_precise:
                precise_a.append(v_item.breakA)
            if v_item.breakB_precise:
                precise_b.append(v_item.breakB)

            informative.append(v_item)

    if not informative:
        return {}

    call_informative = Counter([(i.svtype, i.join_type) for i in informative]).most_common()

    cdef str svtype, jointype
    svtype, jointype = call_informative[0][0]
    dc = make_call(informative, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev)

    return dc


cdef dict one_edge(infile, list u_reads, list v_reads, int clip_length, int insert_size, int insert_stdev,
                   int min_support, int block_edge=1, assemble=True):

    attrs = count_attributes(u_reads, v_reads, min_support)
    if not attrs:
        return {}

    info = call_from_reads(u_reads, v_reads, insert_size, insert_stdev)

    if not info:
        return {}

    info.update(attrs)

    as1 = None
    as2 = None
    if assemble:

        if info["preciseA"]:
            as1 = assembler.base_assemble(u_reads)

        if info["preciseB"]:
            as2 = assembler.base_assemble(v_reads)

    info["linked"] = 0
    info["block_edge"] = block_edge
    info["contig"] = None
    info["contig2"] = None
    rbases = 0
    if as1 is not None and "contig" in as1:
        info["contig"] = as1["contig"]
        rbases += as1["ref_bases"]

    if as2 is not None and "contig" in as2:
        info["contig2"] = as2["contig"]
        rbases += as2["ref_bases"]
    info["ref_bases"] = rbases

    return info


def multi(dict data, bam, int insert_size, int insert_stdev, int clip_length, int min_support):

    # Sometimes partitions are not linked, happens when there is not much support between partitions
    # Then need to decide whether to call from a single partition

    seen = set(data["parts"].keys())

    out_counts = defaultdict(int)  # The number of 'outward' links to other clusters
    for (u, v), d in data["s_between"].items():

        input_reads = fetch_reads(data, d, bam)  # {Node: alignment,..}

        #     echo("main edge")
        #     for nodename in d[u]:
        #         echo(nodename, str(input_reads[nodename]))
        #     for nodename in d[v]:
        #         echo(nodename, str(input_reads[nodename]))
        rd_u = []
        for n in d[u]:
            try:
                rd_u.append(input_reads[n])
            except KeyError:
                echo("Warning missing u read", n)

        out_counts[u] += len(rd_u)

        rd_v = []
        for n in d[v]:
            try:
                rd_v.append(input_reads[n])
            except KeyError:
                echo("Warning missing v read", n)

        out_counts[v] += len(rd_v)

        if len(rd_u) == 0 or len(rd_v) == 0:
            continue

        if u in seen:
            seen.remove(u)
        if v in seen:
            seen.remove(v)
        # if 5371604 in d[u] or 5371604 in d[v]:
        #     echo("u, v", u, v)
        #     echo(one_edge(bam, rd_u, rd_v, clip_length, insert_size, insert_stdev))

        yield one_edge(bam, rd_u, rd_v, clip_length, insert_size, insert_stdev, min_support)

    # Process any unconnected blocks
    if seen:
        for part_key in seen:
            d = data["parts"][part_key]
            if len(d) >= min_support:
                # Call single block, only collect local reads to the block
                rds = {}
                to_collect = {}
                for v in d:
                    try:
                        rds[v] = data["reads"][v]  # May have been collected already
                    except KeyError:
                        to_collect[v] = data["n2n"][v]

                rds.update(graph.get_reads(bam, to_collect))

                if len(rds) < min_support:
                    continue

                yield single(bam, {"reads": rds}, insert_size, insert_stdev, clip_length, min_support, to_assemble=True)

    # Check for events within clustered nodes - happens rarely
    for k, d in data["s_within"].items():

        o_count = out_counts[k]
        i_counts = len(d)
        if o_count > 0 and i_counts > (2*min_support) and i_counts > o_count:

            rds = {}
            to_collect = {}
            for v in d:
                try:
                    rds[v] = data["reads"][v]  # May have been collected already
                except KeyError:
                    to_collect[v] = data["n2n"][v]

            rds.update(graph.get_reads(bam, to_collect))

            if len(rds) < min_support:
                continue

            yield single(bam, {"reads": rds}, insert_size, insert_stdev, clip_length, min_support, to_assemble=True)


def call_from_block_model(bam, data, clip_length, insert_size, insert_stdev, min_support):

    n_parts = len(data["parts"])
    # n_reads = len(data["reads"])
    # if 279 in data["reads"]:
    # echo("reads and parts", n_reads, n_parts)
    # echo(data["reads"])
    # quit()
    if n_parts >= 1:
        # Processed single edges and break apart connected
        for event in multi(data, bam, insert_size, insert_stdev, clip_length, min_support):
        # if 279 in data["reads"]:
        #     echo("called2", event)
            yield event

    elif n_parts == 0: # and min_support == 1:
        # Possible single read only
        yield single(bam, data, insert_size, insert_stdev, clip_length, min_support, to_assemble=1)


cpdef dict get_raw_coverage_information(r, regions, regions_depth, infile):

    # Check if side A in regions
    ar = False  # c_io_funcs.intersecter_int_chrom
    if data_io.intersecter(regions, r["chrA"], r["posA"], r["posA"] + 1):
        ar = True

    if "chrB" not in r:  # Todo Does this happen?
        return None

    br = False
    if data_io.intersecter(regions, r["chrB"], r["posB"], r["posB"] + 1):
        br = True

    # Put non-region first
    kind = None

    if not ar and not br:
        kind = "extra-regional"
        if r["chrA"] == r["chrB"] and r["posA"] > r["posB"]:
            switch = True

        # Skip if regions have been provided; almost always false positives?
        # if regions is not None:
        #     return None

    switch = False
    if (br and not ar) or (not br and ar):
        kind = "hemi-regional"
        if not br and ar:
            switch = True

    if ar and br:

        if r["chrA"] == r["chrB"]:
            rA = list(regions[r["chrA"]].find_overlap(r["posA"], r["posA"] + 1))[0]
            rB = list(regions[r["chrB"]].find_overlap(r["posB"], r["posB"] + 1))[0]

            if rA[0] == rB[0] and rA[1] == rB[1]:
                kind = "intra_regional"
                # Put posA first
                if r["posA"] > r["posB"]:
                    switch = True

            else:
                kind = "inter-regional"
                if r["chrA"] != sorted([r["chrA"], r["chrB"]])[0]:
                    switch = True
        else:
            kind = "inter-regional"

    if switch:
        chrA, posA, cipos95A, contig2 = r["chrA"], r["posA"], r["cipos95A"], r["contig2"]
        r["chrA"] = r["chrB"]
        r["posA"] = r["posB"]
        r["cipos95A"] = r["cipos95B"]
        r["chrB"] = chrA
        r["posB"] = posA
        r["cipos95B"] = cipos95A
        r["contig2"] = r["contig"]
        r["contig"] = contig2

    if kind == "hemi-regional":
        chrom_i = infile.get_tid(r["chrA"])
        if chrom_i in regions_depth:
            reads_10kb = coverage.calculate_coverage(r["posA"] - 10000, r["posA"] + 10000, regions_depth[chrom_i])
            reads_10kb = round(reads_10kb, 3)
        else:
            reads_10kb = 0
    else:
        # Calculate max
        chrom_i = infile.get_tid(r["chrA"])
        if chrom_i in regions_depth:
            reads_10kb_left = coverage.calculate_coverage(r["posA"] - 10000, r["posA"] + 10000, regions_depth[chrom_i])
            reads_10kb_left = round(reads_10kb_left, 3)
        else:
            reads_10kb_left = 0

        chrom_i = infile.get_tid(r["chrB"])
        if chrom_i in regions_depth:
            reads_10kb_right = coverage.calculate_coverage(r["posB"] - 10000, r["posB"] + 10000, regions_depth[chrom_i])
            reads_10kb_right = round(reads_10kb_right, 3)
        else:
            reads_10kb_right = 0

        if reads_10kb_left > reads_10kb_right:
            reads_10kb = reads_10kb_left
        else:
            reads_10kb = reads_10kb_right

    r["kind"] = kind
    r["raw_reads_10kb"] = reads_10kb
    if r["chrA"] != r["chrB"]:
        r["svlen"] = 1000000
    else:
        r["svlen"] = abs(r["posB"] - r["posA"])
    r["su"] = r["pe"] + r["supp"]
    return r
