# cython: language_level=3

from __future__ import absolute_import
from collections import Counter, defaultdict
import itertools
import click
import numpy as np
from . import data_io, assembler, graph, coverage
import pandas as pd
from subprocess import call
import os
import time
import pysam


def echo(*args):
    click.echo(args, err=True)


cdef tuple guess_break_point(read, int node=0):
    # Node id is optional

    # If the read is non-discordant and no clip is present, skip
    cdef int p, t, left, right, exact, flag, rname
    flag = read.flag

    ct = read.cigartuples
    rname = read.rname
    if flag & 2 and not any(i[0] == 4 or i[0] == 5 for i in ct):
        exact = 0
        # Breakpoint position is beyond the end of the last read
        if flag & 16:  # Read on reverse strand guess to the left
            p = read.pos
            t = 5
        else:  # Guess right
            p = read.reference_end
            t = 3

        return t, rname, p, 0, flag, exact, node  #, read.rnext, read.pnext

    # Sometimes a clip may be present, use this as a break point if available
    left = 0
    right = 0
    if ct[0][0] == 4 or ct[0][0] == 5:  # Left soft or hard-clip
        left = ct[0][1]
    if ct[-1][0] == 4 or ct[-1][0] == 5:
        right = ct[-1][0]

    if left > 0 or right > 0:
        exact = 1
        if left > right:
            # prime, rname, pos, soft-clipped, flag, exact
            return 5, rname, read.pos, 1, flag, exact, node  #, read.rnext, read.pnext
        else:
            return 3, rname, read.reference_end, 1, flag, exact, node  #, read.rnext, read.pnext
    else:
        exact = 0
        # Breakpoint position is beyond the end of the last read
        if flag & 16:  # Read on reverse strand guess to the left
            p = read.pos
            t = 5
        else:  # Guess right
            p = read.reference_end
            t = 3
        return t, rname, p, 0, flag, exact, node  #, read.rnext, read.pnext


cpdef dict call_break_points(list c1, list c2):
    """
    Makes a call from a list of break points. Can take a list of break points, or one merged cluster of breaks.
    Outliers are dropped. Breakpoints are clustered into sets
    :param c1: A 5 tuple (3 or 5 join, chromosome, break point position, soft-clipped)
    :param c2: Same as c1
    :return: Info dict containing a summary of the call
    """

    info = {}
    cdef int count = 0
    cdef str x, y, side
    cdef int mean_pos, mean_95, chrom, sc_side
    cdef list grp
    cdef tuple t

    # Fields are: prime, rname, pos, soft-clipped, flag, exact
    for grp in (c1, c2):

        if len(grp) == 0:
            continue  # When no c2 is found

        if count == 0:
            side = "A"
        else:
            side = "B"

        f_strand, r_strand = 0, 0
        for item in grp:
            if item[0] == 3:
                f_strand += 1
            else:
                r_strand += 1
        info["strand" + side] = (f_strand, r_strand)

        soft_clipped = [i for i in grp if i[3]]

        if len(soft_clipped) > len(grp) / 2.:
            grp = soft_clipped  # Often true for small variants

        if len(grp) == 1:
            t = grp[0]
            chrom = t[1]
            sc_side = t[0]
            mean_pos = t[2]
            mean_95 = 0

        else:
            chrom = Counter([i[1] for i in grp]).most_common()[0][0]
            grp = [i for i in grp if i[1] == chrom]

            sc_side = Counter([i[0] for i in grp]).most_common()[0][0]
            bp = [i[2] for i in grp]
            mean_pos = int(np.mean(bp))
            mean_95 = abs(int(np.percentile(bp, [97.5])) - mean_pos)


        info["chr" + side] = chrom
        info["pos" + side] = mean_pos
        info["cipos95" + side] = mean_95
        info["join" + side] = sc_side
        count += 1

    if "chrB" in info and info["chrA"] == info["chrB"]:
        if info["joinA"] == info["joinB"]:
            info["svtype"] = "INV"
            info["join_type"] = f"{info['joinA']}to{info['joinB']}"
        else:
            if info["posA"] <= info["posB"]:
                x = "joinA"
                y = "joinB"
            else:
                x = "joinB"
                y = "joinA"
            if info[x] == 3 and info[y] == 5:
                info["svtype"] = "DEL"
                info["join_type"] = "3to5"
            elif info[x] == 5 and info[y] == 3:
                info["svtype"] = "DUP"
                info["join_type"] = "5to3"

    elif "chrB" in info:
        info["svtype"] = "TRA"
        info["join_type"] =  f"{info['joinA']}to{info['joinB']}"
    else:
        info["svtype"] = "BND"
        if "joinA" in info:
            info["join_type"] = f"{info['joinA']}to?"
        else:
            info["join_type"] = "?to?"

    return info


cdef dict count_attributes(list reads1, list reads2):

    if len(reads1) == 0:
        raise ValueError("No reads in set")

    r = {"pe": 0, "supp": 0, "sc": 0, "DP": [], "DApri": [], "DN": [], "NMpri": [], "NP": 0, "DAsupp": [], "NMsupp": [],
         "maxASsupp": [], "MAPQpri": [], "MAPQsupp": []}
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

            ct = a.cigartuples
            if ct[0][0] == 4 or ct[-1][0] == 4:
                r["sc"] += 1

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


def guess_break_point_from_mate(reads, candidates, max_dist, infile):
    # Candidates is the already processed list from guess_break_point, supplement this with SA information
    # Make a list of possible breaks in the form;
    # (node_name, i), where i = prime, rname, pos, soft-clipped, flag, exact

    count = 0
    for node_name, r in reads.items():

        # If the read is non-discordant and no clip is present, skip
        flag = r.flag
        chrom = r.rname
        rnext = r.rnext
        pnext = r.pnext
        if chrom == rnext and abs(pnext - r.pos) < max_dist:  # Same loci

            if r.has_tag("SA"):  # Parse SA, first alignment is the other read primary line
                sa = r.get_tag("SA").split(",", 3)
                chrom2 = infile.gettid(sa[0])
                pos2 = int(sa[1])
                p1_neg = bool(flag & 16)
                p2_neg = sa[3] == "-"
                same_strand = (p1_neg and p2_neg) or (not p1_neg and not p2_neg)
                if same_strand:
                    if p1_neg:
                        strand = 3
                    else:
                        strand = 5
                else:
                    if p1_neg:
                        strand = 5
                    else:
                        strand = 3

                if chrom2 != chrom or abs(pos2 - r.pos) > (max_dist * 3):

                    candidates.append((node_name, (strand, chrom2, pos2, False, flag, False)))
                    count += 1
        else:
            candidates.append((node_name, (0, rnext, pnext, False, flag, False)))
            count += 1
    if count == 0:
        return None

    return candidates


def one_read(infile, data, insert_size, insert_stdev, clip_length):
    # Same as single (below) but with no assembly step
    return single(infile, data, insert_size, insert_stdev, clip_length, 0, to_assemble=False)


cpdef tuple BFS_local(G, int source, int thresh):

    # Mark all the vertices as not visited
    visited = set([])

    # Create a queue for BFS
    queue = [source]
    nodes_found = set([])
    cdef int u, v

    while queue:
        u = queue.pop(0)

        for v in G.neighbors(u):

            if v not in visited:
                if G.weight(u, v) < thresh:
                    if u not in nodes_found:
                        nodes_found.add(u)
                    if v not in nodes_found:
                        nodes_found.add(v)
                        queue.append(v)

        visited.add(u)

    return nodes_found, visited


cdef class ReadCluster:
    cdef public list breaks
    cdef public int chrom
    # cdef public int chrom_next
    # cdef public int mean_pos_next
    cdef public int mean_pos
    cdef public int count

    def __init__(self, vals):
        self.breaks = vals
        self.chrom = vals[0][2]
        self.mean_pos = int(np.median([i[2] for i in vals]))
        # self.chrom_next = vals[0][7]
        # self.mean_pos_next = int(np.median([i[8] for i in vals]))
        self.count = len(vals)

    cpdef int distance_to(self, int chrom, int pos, int thresh):
        cdef int dis
        if chrom == self.chrom:
            dis = int(abs(self.mean_pos - <float> pos))
            if dis < thresh:
                return dis
        return -1

    cpdef void add_item(self, v):
        self.breaks.append(v)
        self.count += 1


cdef list cluster_by_distance(list bpt, int t, int t2):
    """
    Naively merge breakpoints into clusters based on a distance threshold.
    :param bpt: list of break tuples
    :param t: clustering distance threshold, if this fails use t2
    :param t2: more relaxed clustering distance
    :return: a list of clusters, each cluster is a dict with items of key=(read_name, flag), value=breakpoint info
    """

    clst = []
    current = []
    inexact = []
    # tuple is join-type, rname, pos, soft-clipped, flag, exact, node-name

    cdef int chrom, pos, exact, last_chrom, last_pos, j

    for i in sorted(bpt, key=lambda x: (x[1], x[2])):  # Sort by chrom and pos

        chrom, pos, exact = i[1], i[2], i[5]
        if len(current) == 0 and exact:
            current.append(i)
            continue

        if exact:
            last_chrom, last_pos = current[-1][1], current[-1][2]
            if chrom == last_chrom and abs(last_pos - pos) < t:
                current.append(i)
            else:
                clst.append(current)
                current = [i]
        else:
            inexact.append(i)

    # Last value
    if len(current) > 0:
        clst.append(current)

    if len(clst) == 2:
        return clst

    if len(inexact) == 0:
        if len(clst) > 2:
            # Choose biggest, make call from one side
            return [sorted(clst, key=lambda x: len(x))[-1]]

        else:
            return clst

    # Greedy cluster
    clusters = [ReadCluster(i) for i in clst]

    # Merge inexact, or start new cluster
    cdef int idx, min_d, distance, break_next
    for item in inexact:
        if len(clusters) == 0:
            clusters.append(ReadCluster([item]))
            continue
        else:
            idx = -1
            min_d = 1000000000
            break_next = 0
            for j, cl in enumerate(clusters):

                distance = cl.distance_to(item[1], item[2], thresh=t2)
                if distance!= -1 and distance < min_d:
                    min_d = distance
                    idx = j
                    if break_next:  # Clusters are sorted so only check neighbour cluster
                        break
                    break_next = 1

            if idx != -1:
                clusters[idx].add_item(item)
            else:
                clusters.append(ReadCluster([item]))

    if len(clusters) == 2:
        return [i.breaks for i in clusters]

    # Return largest 2 clusters, doesnt always make sense; cluster 1 and 2 could involve different breaks
    # works well enough most of thr time
    return [i.breaks for i in sorted(clusters, key=lambda x: x.count, reverse=True)][:2]


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
    for temp_name, alignments in tmp.items():
        l_align = list(alignments)
        if len(l_align) > 1:
            pair = guess_informative_pair(l_align)
            # if temp_name == "HWI-D00360:5:H814YADXX:2:2104:14122:33949":
            #     echo(pair[0].to_string())
            #     echo(pair[1].to_string())
            if pair:
                u_reads.append(pair[0])
                v_reads.append(pair[1])
                process_pair(pair, precise_a, precise_b, informative)

    if not informative:
        return {}

    call_informative = Counter([(i["svtype"], i["join_type"]) for i in informative]).most_common()
    svtype, jointype = call_informative[0][0]
    info = make_call(informative, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev)

    if not info:
        return {}

    as1 = None
    as2 = None

    if info["preciseA"]:
        as1 = assembler.base_assemble(u_reads)

    if info["preciseB"]:
        as2 = assembler.base_assemble(v_reads)

    info["linked"] = 0

    info.update(count_attributes(u_reads, v_reads))

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


cdef dict single_old(infile, dict data, int insert_size, int insert_stdev, int clip_length, int min_support,
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

    query_breaks = [guess_break_point(r, n) for n, r in data["reads"].items()]

    if not query_breaks or not len(query_breaks) >= min_support:
        return {}

    if len(query_breaks) == 2:
        u_reads = [data["reads"][query_breaks[0][6]]]
        v_reads = [data["reads"][query_breaks[1][6]]]

        info = one_edge(infile, u_reads, v_reads, clip_length, insert_size, insert_stdev,
                        block_edge=0, assemble=to_assemble)
        return info

    clst = cluster_by_distance(query_breaks, t=25, t2=insert_stdev + (2*insert_stdev))

    if len(clst) == 2:
        # echo("0", clst[0])
        # echo("1", clst[1])

        c1, c2 = clst

        u_reads = [data["reads"][t[6]] for t in c1]
        v_reads = [data["reads"][t[6]] for t in c2]

        info = one_edge(infile, u_reads, v_reads, clip_length, insert_size, insert_stdev,
                        block_edge=0, assemble=to_assemble)

        # if 279 in data["reads"]:
        #     echo("single cll", len(u_reads), len(v_reads))
        #     echo(info)
        #
        #     echo(call_from_reads(u_reads, v_reads, insert_size, insert_stdev))
        return info

    elif len(clst) == 1:
        c1 = clst[0]  # Assume one cluster
        pass

    return {}
#
# cdef tuple get_tuple(dict j):
#     # Get breakpoint info from contig. Need (3 or 5 join, read names, chromosome, break point position, soft-clipped)
#     return (5 if j["left_clips"] > j["right_clips"] else 3,
#             j["bamrname"],
#             j["ref_start"] if j["left_clips"] > j["right_clips"] else j["ref_end"],
#             j["contig"][0].islower() or j["contig"][-1].islower())


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


cdef void two_primary(dict d):

    if d["posA"] < d["posB"] or (d["posA"] == d["posB"] and d["endA"] < d["endB"]):

        if d["strandA"] == 3 and d["strandB"] == 5:  # DEL type
            d["breakA"] = d["endA"]
            if d["right_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["posB"]
            if d["left_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DEL"
            d["join_type"] = "3to5"

        elif d["strandA"] == 5 and d["strandB"] == 3:  # DUP type
            d["breakA"] = d["posA"]
            if d["left_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["endB"]
            if d["right_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DUP"
            d["join_type"] = "5to3"

        elif d["strandA"] == d["strandB"]:  # INV type
            # Break to left
            if d["strandA"] == 5:

                if not (d["left_clipA"] or d["left_clipB"]) and (d["right_clipA"] or d["right_clipB"]):
                    d["breakA"] = d["endA"]
                    if d["right_clipA"]:
                        d["breakA_precise"] = 1

                    d["breakB"] = d["endB"]
                    if d["right_clipB"]:
                        d["breakB_precise"] = 1

                    d["svtype"] = "INV"
                    d["join_type"] = "3to3"

                else:
                    d["breakA"] = d["posA"]
                    if d["left_clipA"]:
                        d["breakA_precise"] = 1

                    d["breakB"] = d["posB"]
                    if d["right_clipB"]:
                        d["breakB_precise"] = 1

                    d["svtype"] = "INV"
                    d["join_type"] = "5to5"

            elif d["strandA"] == 3:

                if not (d["right_clipA"] or d["right_clipB"]) and (d["left_clipA"] or d["left_clipB"]):
                    d["breakA"] = d["posA"]
                    if d["right_clipA"]:
                        d["breakA_precise"] = 1

                    d["breakB"] = d["posB"]
                    if d["right_clipB"]:
                        d["breakB_precise"] = 1

                    d["svtype"] = "INV"
                    d["join_type"] = "5to5"

                else:

                    d["breakA"] = d["endA"]
                    if d["left_clipA"]:
                        d["breakA_precise"] = 1

                    d["breakB"] = d["endB"]
                    if d["right_clipB"]:
                        d["breakB_precise"] = 1

                    d["svtype"] = "INV"
                    d["join_type"] = "3to3"

            # if not d["right_clipA"] and not d["right_clipB"]:
            #     d["breakA"] = d["posA"]
            #     if d["left_clipA"]:
            #         d["breakA_precise"] = 1
            #
            #     d["breakB"] = d["posB"]
            #     if d["left_clipB"]:
            #         d["breakB_precise"] = 1
            #
            #     d["svtype"] = "INV"
            #     d["join_type"] = "5to5"

            # elif d["strandA"] == 5:
            #     d["breakA"] = d["posA"]
            #     if d["left_clipA"]:
            #         d["breakA_precise"] = 1
            #
            #     d["breakB"] = d["posB"]
            #     if d["right_clipB"]:
            #         d["breakB_precise"] = 1
            #
            #     d["svtype"] = "INV"
            #     d["join_type"] = "5to5"

            # Break to right
            # elif not d["left_clipA"] and not d["left_clipB"]:
            #     d["breakA"] = d["endA"]
            #     if d["right_clipA"]:
            #         d["breakA_precise"] = 1
            #
            #     d["breakB"] = d["endB"]
            #     if d["right_clipB"]:
            #         d["breakB_precise"] = 1
            #
            #     d["svtype"] = "INV"
            #     d["join_type"] = "3to3"

            # elif d["strandA"] == 3:
            #     d["breakA"] = d["endA"]
            #     if d["left_clipA"]:
            #         d["breakA_precise"] = 1
            #
            #     d["breakB"] = d["endB"]
            #     if d["right_clipB"]:
            #         d["breakB_precise"] = 1
            #
            #     d["svtype"] = "INV"
            #     d["join_type"] = "3to3"

            # Break INV:DUP
            else:
            # elif d["left_clipA"] and d["right_clipB"]:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV:DUP"
                d["join_type"] = "3to5"

    else:  # B < A

        if d["strandA"] == 5 and d["strandB"] == 3:  # DEL type
            d["breakA"] = d["posA"]
            if d["left_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["endB"]
            if d["right_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DEL"
            d["join_type"] = "3to5"

        elif d["strandB"] == 5 and d["strandA"] == 3:  # DUP type
            d["breakA"] = d["endA"]
            if d["right_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["posB"]
            if d["left_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DUP"
            d["join_type"] = "5to3"

        # else:  # INV type
            # Break to left
            # if not d["right_clipA"] and not d["right_clipB"]:
            #     d["breakA"] = d["posA"]
            #     if d["left_clipA"]:
            #         d["breakA_precise"] = 1
            #     d["breakB"] = d["posB"]
            #     if d["left_clipB"]:
            #         d["breakB_precise"] = 1
            #     d["svtype"] = "INV"
            #     d["join_type"] = "5to5"
            #
            # # Break to the right
            # elif not d["left_clipA"] and not d["left_clipB"]:
            #     d["breakA"] = d["endA"]
            #     if d["right_clipA"]:
            #         d["breakA_precise"] = 1
            #     d["breakB"] = d["endB"]
            #     if d["right_clipB"]:
            #         d["breakB_precise"] = 1
            #     d["svtype"] = "INV"
            #     d["join_type"] = "3to3"
            #
            # elif d["left_clipA"] and d["right_clipB"]:
            #     d["breakA"] = d["posA"]
            #     d["breakA_precise"] = 1
            #     d["breakB"] = d["endB"]
            #     d["breakB_precise"] = 1
            #     d["svtype"] = "INV:DUP"
            #     d["join_type"] = "5to3"
            #
            # else:
            # # elif d["right_clipA"] and d["left_clipB"]:
            #     d["breakA"] = d["endA"]
            #     d["breakA_precise"] = 1
            #     d["breakB"] = d["posB"]
            #     d["breakB_precise"] = 1
            #     d["svtype"] = "INV:DUP"
            #     d["join_type"] = "3to5"
            # INV type
        elif d["strandA"] == 5:

            if not (d["left_clipA"] or d["left_clipB"]) and (d["right_clipA"] or d["right_clipB"]):
                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            else:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["posB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV"
                d["join_type"] = "5to5"

        elif d["strandA"] == 3:

            if not (d["right_clipA"] or d["right_clipB"]) and (d["left_clipA"] or d["left_clipB"]):
                d["breakA"] = d["posA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["posB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV"
                d["join_type"] = "5to5"

            else:

                d["breakA"] = d["endA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV"
                d["join_type"] = "3to3"


cdef void same_read(dict d):

    # A is first
    if d["posA"] < d["posB"] or (d["posA"] == d["posB"] and d["endA"] <= d["endB"]):

        if d["strandA"] == d["strandB"]:  # Same for 3 and 5 strand reads

            if is_overlapping(d["posA"], d["endA"], d["posB"], d["endB"]):  # Nested DUP
                if d["left_clipA"]:
                    d["breakA"] = d["posA"]
                    d["breakA_precise"] = 1
                elif d["right_clipA"]:
                    d["breakA"] = d["endA"]
                    d["breakA_precise"] = 1
                else:
                    d["breakA"] = d["endA"]

                if d["left_clipB"]:
                    d["breakB"] = d["posB"]
                    d["breakB_precise"] = 1
                elif d["right_clipB"]:
                    d["breakB"] = d["endB"]
                    d["breakB_precise"] = 1
                else:
                    d["breakB"] = d["endB"]
                d["svtype"] = "DUP"
                d["join_type"] = "5to3"

            elif not d["left_clipA"] and not d["right_clipB"]:
                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["posB"]
                if d["left_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "DEL"
                d["join_type"] = "3to5"

            elif not d["right_clipA"] and not d["left_clipB"]:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "DUP"
                d["join_type"] = "5to3"

            elif not d["left_clipA"] and not d["left_clipB"]:
                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            else:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["posB"]
                if d["left_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "5to5"

        elif d["strandA"] == 5 and d["strandB"] == 3:

            # Break right
            if not d["left_clipA"] and not d["left_clipB"]:
                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            # Break left
            elif not d["right_clipA"] and not d["right_clipB"]:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1

                d["breakB"] = d["posB"]
                if d["left_clipB"]:
                    d["breakB_precise"] = 1

                d["svtype"] = "INV"
                d["join_type"] = "5to5"

            elif is_overlapping(d["posA"], d["endA"], d["posB"], d["endB"]):  # Inverted duplication
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV:DUP"
                d["join_type"] = "5to3"

            elif d["right_clipA"] and d["left_clipB"]:
                d["breakA"] = d["endA"]
                d["breakB"] = d["posB"]
                d["svtype"] = "INV"
                d["join_type"] = "3to5"

            else:
            # elif d["left_clipA"] and d["right_clipB"]:
                d["breakA"] = d["posA"]
                d["breakB"] = d["endB"]
                d["svtype"] = "INV"
                d["join_type"] = "5to3"

        else:  # INV type
            # Break left
            if d["left_clipA"] and d["left_clipB"]:
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["breakA"] = d["posA"]
                d["breakB"] = d["posB"]
                d["svtype"] = "INV"
                d["join_type"] = "5to5"
            # Break right
            elif d["right_clipA"] and d["right_clipB"]:
                d["breakA"] = d["endA"]
                d["breakB"] = d["endB"]
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            else:  # Guess using pos only
                d["breakA"] = d["posA"]
                d["breakB"] = d["posB"]
                d["svtype"] = "INV"
                if d["strandA"] == 5:
                    d["join_type"] = "5to5"
                else:
                    d["join_type"] = "3to3"

    # B is first
    else:
        if d["strandA"] == d["strandB"]:

            if is_overlapping(d["posA"], d["endA"], d["posB"], d["endB"]):  # Nested DUP
                if d["left_clipA"]:
                    d["breakA"] = d["posA"]
                    d["breakA_precise"] = 1
                elif d["right_clipA"]:
                    d["breakA"] = d["endA"]
                    d["breakA_precise"] = 1
                else:
                    d["breakA"] = d["endA"]

                if d["left_clipB"]:
                    d["breakB"] = d["posB"]
                    d["breakB_precise"] = 1
                elif d["right_clipB"]:
                    d["breakB"] = d["endB"]
                    d["breakB_precise"] = 1
                else:
                    d["breakB"] = d["endB"]
                d["svtype"] = "DUP"
                d["join_type"] = "5to3"

            elif not d["left_clipB"] and not d["right_clipA"]:
                d["breakA"] = d["posA"]
                d["breakB"] = d["endB"]
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["svtype"] = "DEL"
                d["join_type"] = "3to5"

            elif not d["right_clipB"] and not d["left_clipA"]:
                d["breakA"] = d["endA"]
                d["breakB"] = d["posB"]
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["svtype"] = "DUP"
                d["join_type"] = "5to3"

            else:
                d["breakA"] = d["posA"]
                d["breakB"] = d["posB"]
                d["svtype"] = "BND"
                d["join_type"] = f"{d['strandA']}to{d['strandB']}"

        else:  # INV type
            # Break left
            if d["left_clipA"] and d["left_clipB"]:
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["breakA"] = d["posA"]
                d["breakB"] = d["posB"]
                d["svtype"] = "INV"
                d["join_type"] = "5to5"

            # Break right
            elif d["right_clipA"] and d["right_clipB"]:
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["breakA"] = d["endA"]
                d["breakB"] = d["endB"]
                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            else:  # Guess using pos only
                d["breakA"] = d["posA"]
                d["breakB"] = d["posB"]
                d["svtype"] = "INV"
                if d["strandA"] == 5:
                    d["join_type"] = "5to5"
                else:
                    d["join_type"] = "3to3"


cdef void different_read(dict d):

    if d["posA"] < d["posB"] or (d["posA"] == d["posB"] and d["endA"] <= d["endB"]):  # A is first

        if d["strandA"] == 3 and d["strandB"] == 5:  # DEL type
            d["breakA"] = d["endA"]
            if d["right_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["posB"]
            if d["left_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DEL"
            d["join_type"] = "3to5"

        elif d["strandA"] == 5 and d["strandB"] == 3:  # DUP type
            d["breakA"] = d["posA"]
            if d["left_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["endB"]
            if d["right_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DUP"
            d["join_type"] = "5to3"

        elif d["strandA"] == d["strandB"]:
            if is_overlapping(d["posA"], d["endA"], d["posB"], d["endB"]):  # Nested

                if d["strandA"] == 3:  # Both forward strand
                    d["breakA"] = d["endA"]
                    if d["right_clipA"]:
                        d["breakA_precise"] = 1

                    d["breakB"] = d["endB"]
                    if d["right_clipB"]:
                        d["breakB_precise"] = 1
                    d["svtype"] = "INV"
                    d["join_type"] = "3to3"

                else:  # Both forward strand
                    d["breakA"] = d["posA"]
                    if d["left_clipA"]:
                        d["breakA_precise"] = 1

                    d["breakB"] = d["posB"]
                    if d["left_clipB"]:
                        d["breakB_precise"] = 1
                    d["svtype"] = "INV"
                    d["join_type"] = "5to5"

            elif not d["right_clipA"] and not d["right_clipB"]:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["posB"]
                if d["left_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "5to5"

            elif not d["left_clipA"] and not d["left_clipB"]:
                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            elif d["right_clipA"] and d["left_clipB"]:
                d["breakA"] = d["endA"]
                d["breakA_precise"] = 1
                d["breakB"] = d["posB"]
                d["breakB_precise"] = 1
                d["svtype"] = "INV:DUP"
                d["join_type"] = "5to3"

            else:
            # elif d["left_clipA"] and d["right_clipB"]:
                d["breakA"] = d["posA"]
                d["breakA_precise"] = 1
                d["breakB"] = d["endB"]
                d["breakB_precise"] = 1
                d["svtype"] = "INV:DUP"
                d["join_type"] = "3to5"

    else:  # B is first; B <= A

        if d["strandA"] == 5 and d["strandB"] == 3:  # DEL type
            d["breakA"] = d["posA"]
            if d["left_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["endB"]
            if d["right_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DEL"
            d["join_type"] = "3to5"

        elif d["strandA"] == 3 and d["strandB"] == 5:  # DUP type
            d["breakA"] = d["endA"]
            if d["right_clipA"]:
                d["breakA_precise"] = 1

            d["breakB"] = d["posB"]
            if d["left_clipB"]:
                d["breakB_precise"] = 1

            d["svtype"] = "DUP"
            d["join_type"] = "5to3"

        elif d["strandA"] == d["strandB"]:  # INV type

            if is_overlapping(d["posA"], d["endA"], d["posB"], d["endB"]):  # Nested DUP
                if d["left_clipA"]:
                    d["breakA"] = d["posA"]
                    d["breakA_precise"] = 1
                elif d["right_clipA"]:
                    d["breakA"] = d["endA"]
                    d["breakA_precise"] = 1
                else:
                    d["breakA"] = d["endA"]

                if d["left_clipB"]:
                    d["breakB"] = d["posB"]
                    d["breakB_precise"] = 1
                elif d["right_clipB"]:
                    d["breakB"] = d["endB"]
                    d["breakB_precise"] = 1
                else:
                    d["breakB"] = d["endB"]
                d["svtype"] = "DUP"
                d["join_type"] = "5to3"

            elif d["right_clipA"] and d["left_clipB"]:
                d["breakA"] = d["endA"]
                d["breakB"] = d["posB"]
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["svtype"] = "DEL"
                d["join_type"] = "3to5"

            elif d["left_clipA"] and d["right_clipB"]:
                d["breakA"] = d["posA"]
                d["breakB"] = d["endB"]
                d["breakA_precise"] = 1
                d["breakB_precise"] = 1
                d["svtype"] = "DUP"
                d["join_type"] = "5to3"

            elif not d["left_clipB"] and not d["left_clipA"]:

                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1

                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            elif not d["right_clipB"] and not d["right_clipA"]:

                d["breakB"] = d["posB"]
                if d["left_clipB"]:
                    d["breakB_precise"] = 1

                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "5to5"

            elif d["strandA"] == 3:
                d["breakA"] = d["endA"]
                if d["right_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["endB"]
                if d["right_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "3to3"

            else:
                d["breakA"] = d["posA"]
                if d["left_clipA"]:
                    d["breakA_precise"] = 1
                d["breakB"] = d["posB"]
                if d["left_clipB"]:
                    d["breakB_precise"] = 1
                d["svtype"] = "INV"
                d["join_type"] = "5to5"


cdef void translocation(dict d):

    d["svtype"] = "TRA"

    if d["left_clipA"]:
        d["breakA"] = d["posA"]
        d["breakA_precise"] = 1
    elif d["right_clipA"]:
        d["breakA"] = d["endA"]
        d["breakA_precise"] = 1
    else:
        d["breakA"] = d["posA"]
        d["breakA_precise"] = 1

    if d["left_clipB"]:
        d["breakB"] = d["posB"]
        d["breakB_precise"] = 1
    elif d["right_clipB"]:
        d["breakB"] = d["endB"]
        d["breakB_precise"] = 1
    else:
        d["breakB"] = d["posB"]
        d["breakB_precise"] = 1
    d["join_type"] = f"{d['strandA']}to{d['strandB']}"


cdef int is_overlapping(int x1, int x2, int y1, int y2) nogil:
    return int(max(x1, y1) <= min(x2, y2))


cdef void classify_d(dict d):

    d["breakA_precise"] = 0
    d["breakB_precise"] = 0

    if d["chrA"] != d["chrB"]:
        translocation(d)

    else:  # Intra-chromosomal
        # Find join type first, different for pri-sup, pri-pri
        if d["priA"] and d["priB"]:  # Both primary
            two_primary(d)

        else:  # One is a supplementary
            if d["rA"] == d["rB"]:
                same_read(d)
            else:
                different_read(d)


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
    positionsA = [i["breakA"] for i in informative]
    positionsB = [i["breakB"] for i in informative]
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

    return {"svtype": svtype, "join_type": jointype, "chrA": informative[0]["chrA"], "chrB": informative[0]["chrB"],
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


cdef void process_pair(tuple pair, list precise_a, list precise_b, list informative):
    a, b = pair
    a_ct = a.cigartuples
    b_ct = b.cigartuples
    # Soft-clips for the chosen pair, plus template start of alignment
    left_clipA, right_clipA, left_clipB, right_clipB = mask_soft_clips(a.flag, b.flag, a_ct, b_ct)

    d = {"chrA": a.rname,
        "chrB": b.rname,
        "priA": not a.flag & 3328,
        "priB": not b.flag & 3328,
        "rA": 1 if a.flag & 64 else 2,
        "rB": 1 if b.flag & 64 else 2,
        "posA": a.pos,
        "endA": a.reference_end,
        "posB": b.pos,
        "endB": b.reference_end,
        "strandA": 3 if a.flag & 16 == 0 else 5,
        "strandB": 3 if b.flag & 16 == 0 else 5,
        "left_clipA": left_clipA,
        "right_clipA": right_clipA,
        "left_clipB": left_clipB,
        "right_clipB": right_clipB,
        "qname": a.qname[-5:],
        "breakA_precise": 0,
        "breakB_precise": 0
        }

    # If both sides clipped choose longest soft-clip (works ok)
    if d["left_clipA"] and d["right_clipA"]:
        if a_ct[0][1] >= a_ct[-1][1]:
            d["right_clipA"] = 0
        else:
            d["left_clipA"] = 0

    if d["left_clipB"] and d["right_clipB"]:
        if b_ct[0][1] >= b_ct[-1][1]:
            d["right_clipB"] = 0
        else:
            d["left_clipB"] = 0

    classify_d(d)

    if d["breakA_precise"]:
        precise_a.append(d["breakA"])
    if d["breakB_precise"]:
        precise_b.append(d["breakB"])

    informative.append(d)


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
    for qname in [k for k in grp_u if k in grp_v]:  # Qname fount on both sides

        u = grp_u[qname]
        v = grp_v[qname]

        if len(u) == 1 and len(v) == 1:
            pair = (u[0], v[0])
        else:
            pair = informative_pair(u, v)

        if pair:
            process_pair(pair, precise_a, precise_b, informative)


    # echo("informative", informative)
    if not informative:
        return {}

    call_informative = Counter([(i["svtype"], i["join_type"]) for i in informative]).most_common()
    svtype, jointype = call_informative[0][0]
    dc = make_call(informative, precise_a, precise_b, svtype, jointype, insert_size, insert_stdev)

    return dc


cdef dict one_edge(infile, list u_reads, list v_reads, int clip_length, int insert_size, int insert_stdev,
                   int block_edge=1, assemble=True):

    info = call_from_reads(u_reads, v_reads, insert_size, insert_stdev)

    if not info:
        return {}

    as1 = None
    as2 = None
    if assemble:

        if info["preciseA"]:
            as1 = assembler.base_assemble(u_reads)

        if info["preciseB"]:
            as2 = assembler.base_assemble(v_reads)

    # if info["posB"] == 46698003:
    #     echo(info)
    #
    #     echo("--")
    #     echo(as1)
    #     echo(as2)

    # if as1 is None or len(as1) == 0:
    #     guess_a = [guess_break_point(r) for r in u_reads]
    # else:
    #     guess_a = [get_tuple(as1)]  # Tuple of breakpoint information

    # if not guess_a:
    #     return {}

    # if as2 is None or len(as2) == 0:
    #     guess_b = [guess_break_point(r) for r in v_reads]
    # else:
    #     guess_b = [get_tuple(as2)]
    #
    # echo()


    # weight_call(u_reads, infile)
    # weight_call(v_reads, infile)


    # if not guess_b:
    #     return {}

    # info = call_break_points(guess_a, guess_b)

    info["linked"] = 0

    info.update(count_attributes(u_reads, v_reads))

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

        yield one_edge(bam, rd_u, rd_v, clip_length, insert_size, insert_stdev)

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
        # echo(o_count > 0, i_counts > (2*min_support), i_counts > (3*o_count), i_counts, o_count)
        if o_count > 0 and i_counts > (2*min_support) and i_counts > (o_count):

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

            # if 279 in d:

            #     echo(len(rds))
            #     echo("here", k)
            yield single(bam, {"reads": rds}, insert_size, insert_stdev, clip_length, min_support, to_assemble=True)



def call_from_block_model(bam, data, clip_length, insert_size, insert_stdev, min_support):

    n_parts = len(data["parts"])
    n_reads = len(data["reads"])
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

    # elif n_parts == -1:
        # Single isolated node, multiple reads

        # d = data["parts"][next(iter(data["parts"]))]  # Get first key

        # rds = {}
        # to_collect = {}
        # for v in d:
        #     if v not in data["reads"]:
            # try:
            #     rds[v] = data["reads"][v]  # May have been collected already

            # except KeyError:
            #     to_collect[v] = data["n2n"][v]

        # rds.update(cy_graph.get_reads(bam, to_collect))
        # data["reads"] = rds
        # data["reads"].update(graph.get_reads(bam, to_collect))

        # if 67 in data["parts"][0]:
        #     echo("single", )
        #     echo(single(bam, data, insert_size, insert_stdev, clip_length, min_support))

        # yield single(bam, data, insert_size, insert_stdev, clip_length, min_support)


    elif n_parts == 0: # and min_support == 1:
        # Possible single read only
        yield single(bam, data, insert_size, insert_stdev, clip_length, 0, to_assemble=1)
    #
    # else:
    #     pass

def get_raw_coverage_information(r, regions, regions_depth, infile):

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
        # Skip if regions have been provided; almost always false positives
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

    return r


def calculate_prob_from_model(all_rows, models):

    if len(all_rows) == 0:
        return

    df = pd.DataFrame.from_records(all_rows).sort_values(["kind", "chrA", "posA"])

    features = ['cipos95A', 'cipos95B', 'DP', 'DApri', 'DN', 'NMpri', 'NP', 'DAsupp', 'NMsupp', 'maxASsupp',
                'contig1_exists', 'both_contigs_exist', 'contig2_exists', 'pe', 'supp', 'sc', 'block_edge', 'MAPQpri',
                'MAPQsupp', 'raw_reads_10kb', 'gc', 'neigh', 'rep', 'ref_bases']
    if not models:
        df["Prob"] = [1] * len(df)  # Nothing to be done
        return df

    df["contig1_exists"] = [1 if str(c) != "None" else 0 for c in df["contig"]]
    df["contig2_exists"] = [1 if str(c) != "None" else 0 for c in df["contig2"]]
    df["both_contigs_exist"] = [1 if i == 1 and j == 1 else 0 for i, j in zip(df["contig1_exists"],
                                                                              df["contig2_exists"])]

    prob = []
    for idx, grp in df.groupby("kind"):

        X = grp[features].astype(float)
        if idx in models:
            clf = models[idx]
            probs = clf.predict_proba(X)
            prob_true = 1 - probs[:, 0]
            prob += list(prob_true)

        else:
            # Choose highest probability out of trained models
            pp = []
            for k in models.keys():
                pp.append(1 - models[k].predict_proba(X)[:, 0])
            a = np.array(pp)
            max_p = np.max(a, axis=0)
            prob += list(max_p)

    df["Prob"] = prob
    return df.sort_values(["kind", "Prob"], ascending=[True, False])
