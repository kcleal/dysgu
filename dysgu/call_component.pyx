# cython: language_level=3

from __future__ import absolute_import
from collections import Counter, defaultdict
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

        grp2 = [i for i in grp if i[3]]
        if len(grp2) > 0:
            grp = grp2

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

        if count == 0:
            side = "A"
        else:
            side = "B"
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
        r["maxASsupp"] = max(r["maxASsupp"])
    else:
        r["maxASsupp"] = 0

    return r


# def pre_process_breakpoints(vals):
#     # If number chroms > 2, reduce
#     if not vals:
#         return []
#     chroms = Counter([i[2] for i in vals if len(i) > 0])
#     if len(chroms) == 0:
#         return []
#
#     if len(chroms) > 2:
#         c = [i[0] for i in sorted(chroms.items(), key=lambda x: x[1], reverse=True)][:2]
#         chroms = {k: v for k, v in chroms.items() if k in c}
#
#     # If minimum chrom counts < 0.1, reduce. Drop low coverage chromosomes
#     if len(chroms) == 2:
#         ci = list(chroms.items())
#         total = float(sum(chroms.values()))
#         if ci[0][1] / total < 0.05:  # Todo add parameter to list
#             del chroms[ci[0][0]]  # Keep item 1
#         elif ci[1][1] / total < 0.05:
#             del chroms[ci[1][0]]  # Keep item 0
#
#     return [v for v in vals if v[2] in chroms]


# def breaks_from_one_side(reads):
#     return [guess_break_point(r) for r in reads]
#     guessed = [guess_break_point(r) for r in reads]
#     return pre_process_breakpoints(guessed)


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


# cdef list get_partitions(G, int thresh, list nodes, int thresh_large):
#
#     seen = set([])  # Dont use Py_IntSet, unions are not quick at the moment
#
#     cdef int u, v, single, i, weight, node
#     cdef float distance
#     parts = []
#     secondary = []
#
#     for u in nodes:
#         if u in seen:
#             continue
#         single = 1
#         for v in G.neighbors(u):
#             if v in seen:
#                 continue
#             weight = G.weight(u, v)
#             if weight < thresh:  # 2 or 3
#                 found, visited_local = BFS_local(G, u, thresh)
#                 seen |= visited_local
#                 parts.append(found)
#                 single = 0
#             elif weight < thresh_large:
#                 secondary.append(u)
#                 single = 0
#
#         if single:
#             parts.append({u})
#
#         seen.add(u)
#
#     if len(parts) == 2:
#         return parts
#
#     if secondary:
#         part_labels = {node: i for i, p in enumerate(parts) for node in p}
#         # echo("labels", part_labels)
#         # echo("parts", parts)
#         for u in secondary:
#             # Add to closest partition, closest neighbor may not be in partition yet
#             closest_vs = sorted([(G.weight(u, v), v) for v in G.neighbors(u)])
#
#             for distance, v in closest_vs:
#                 if v in part_labels:
#                     parts[part_labels[v]].add(u)
#                     break
#
#     return sorted(parts, key=lambda x: len(x), reverse=True)


# cdef list partition_reads(list bpt, int close_thresh, int far_thresh):
#     # echo(bpt)
#     bpt = sorted(bpt, key=lambda x: (x[1], x[2]))
#     partitions = []
#     G = nk.Graph(len(bpt), weighted=True)
#     # tuple is join-type, rname, pos, soft-clipped, flag, exact, node-name
#     cdef int i, dis
#     cdef tuple t1, t2
#
#     for i in range(len(bpt) - 1):
#
#         t1 = bpt[i]
#         t2 = bpt[i+1]
#
#         if t1[1] != t2[1]:
#             continue
#
#         dis = abs(t1[2] - t2[2])
#
#         G.addEdge(i, i+1, w=dis)
#
#     parts = get_partitions(G, close_thresh, list(range(len(bpt))), far_thresh)
#     # longest = float(len(parts[0]))
#     # Drop low coverage partitions
#     # parts = [p for p in parts if len(p) / longest > 0.2]
#
#     return [[bpt[i] for i in p] for p in parts[:2]]


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


cdef dict single(infile, dict data, int insert_size, int insert_stdev, int clip_length, int min_support,
                 int to_assemble=1):
    # Make sure at least one read is worth calling
    cdef int min_distance = insert_size + (2*insert_stdev)
    cdef int n_templates = len(set([i.qname for i in data["reads"].values()]))
    # if n_templates < min_support:
    #     return {}

    if n_templates == 1:
        if not any((not i.flag & 2) or (i.rname != i.rnext) or (abs(i.tlen) > min_distance) for i in data["reads"].values()):
            return {}

    # Single's can sometimes be seperated into two groups --> send one_edge
    # Otherwise, infer the other breakpoint from a single group

    query_breaks = [guess_break_point(r, n) for n, r in data["reads"].items()]

    # if 179 in data["reads"]:
    #     echo("q breaks", len(query_breaks), query_breaks)
    #     quit()
    if not query_breaks or not len(query_breaks) >= min_support:
        return {}

    if len(query_breaks) == 2:
        u_reads = [data["reads"][query_breaks[0][6]]]
        v_reads = [data["reads"][query_breaks[1][6]]]

        info = one_edge(u_reads, v_reads, clip_length, block_edge=0, assemble=to_assemble)
        return info

    clst = cluster_by_distance(query_breaks, t=25, t2=insert_stdev + (2*insert_stdev))
    # if 179 in data["reads"]:
    #     echo(clst)
    if len(clst) == 2:
        # echo(clst)
        c1, c2 = clst
        # echo(c1)
        # echo(c2)
        # echo()
        u_reads = [data["reads"][t[6]] for t in c1]
        v_reads = [data["reads"][t[6]] for t in c2]

        info = one_edge(u_reads, v_reads, clip_length, block_edge=0, assemble=to_assemble)

        # if 41 in data["reads"]:
        #     echo(info)
        return info
        # quit()

    elif len(clst) == 1:
        c1 = clst[0]  # Assume one cluster
        pass

    return {}
    # dict_a, dict_b = separate_mixed(list(breaks_dict.items()), thresh=insert_size + insert_stdev)

    # If an exact breakpoint was found in both, assume they are well seperated and send one_edge
    # found = 0
    # if len(dict_b) > 0:
    #     dav = list(dict_a.values())
    #     dbv = list(dict_b.values())
    #     if any(v[5] for v in dav) and any(v[5] for v in dbv):
    #         found = 1
    #     elif dav[0][1] != dbv[0][1]:  # Different chroms
    #         found = 1
    #     else:  # Make sure regions are well seperated
    #         pa = np.mean([i[2] for i in dav])
    #         pb = np.mean([i[2] for i in dbv])
    #         if abs(pb - pa) > (insert_size + (2*insert_stdev)):
    #             found = 1
    #
    # if found:
    #
    #     u_reads = [data["reads"][n] for n in dict_a]
    #     v_reads = [data["reads"][n] for n in dict_b]
    #
    #     info = one_edge(u_reads, v_reads, clip_length, block_edge=0, assemble=to_assemble)
    #     return info

    # else:
    #
    #     # Infer other locus from mate information
    #     candidates = list(breaks_dict.items())
    #     max_dist = insert_size + (3*insert_stdev)
    #     candidates = guess_break_point_from_mate(data["reads"], candidates, max_dist, infile)
    #
    #     if candidates is not None:
    #
    #         dict_a, dict_b = sorted(separate_mixed(candidates, thresh=max_dist), key=lambda x: len(x), reverse=True)
    #
    #
    #         u_reads = [data["reads"][n] for n in dict_a]
    #         as1 = cy_assembler.base_assemble(u_reads)
    #         if as1:
    #             guess_a = [get_tuple(as1)]
    #         else:
    #             guess_a = [guess_break_point(r) for r in dict_a.values()] #  pre_process_breakpoints(dict_a.values())
    #         guess_b = [guess_break_point(r) for r in dict_b.values()]  # pre_process_breakpoints(dict_b.values())
    #         if not guess_b:
    #             return None
    #
    #         info = call_break_points(guess_a, guess_b)
    #
    #         info["linked"] = 0
    #
    #         info.update(count_attributes(u_reads, []))
    #
    #         info["block_edge"] = -1
    #         info["contig"] = None
    #         info["contig2"] = None
    #
    #         if as1 is not None and "contig" in as1:
    #             info["contig"] = as1["contig"]
    #
    #         return info

cdef tuple get_tuple(dict j):
    # Get breakpoint info from contig. Need (3 or 5 join, read names, chromosome, break point position, soft-clipped)
    return (5 if j["left_clips"] > j["right_clips"] else 3,
            j["bamrname"],
            j["ref_start"] if j["left_clips"] > j["right_clips"] else j["ref_end"],
            j["contig"][0].islower() or j["contig"][-1].islower())


cdef dict one_edge(list u_reads, list v_reads, int clip_length, int block_edge=1, assemble=True):

    if assemble:
        as1 = assembler.base_assemble(u_reads)
        as2 = assembler.base_assemble(v_reads)
    else:
        as1 = None
        as2 = None

    if as1 is None or len(as1) == 0:
        guess_a = [guess_break_point(r) for r in u_reads]
    else:
        guess_a = [get_tuple(as1)]  # Tuple of breakpoint information

    if not guess_a:
        return {}

    if as2 is None or len(as2) == 0:
        guess_b = [guess_break_point(r) for r in v_reads]
    else:
        guess_b = [get_tuple(as2)]
    if not guess_b:
        return {}

    info = call_break_points(guess_a, guess_b)

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

    out_counts = defaultdict(int)
    for (u, v), d in data["s_between"].items():

        input_reads = fetch_reads(data, d, bam)  # {Node: alignment,..}

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

        yield one_edge(rd_u, rd_v, clip_length)

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

    # Check for events within nodes - happens rarely
    for k, d in data["s_within"].items():

        o_count = out_counts[k]
        i_counts = len(d)

        if o_count > 0 and i_counts > (2*min_support) and i_counts > (3*o_count):

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

    # if isinstance(bam, tuple):
    #     bam = pysam.AlignmentFile(bam[0], bam[1])  # Re open
    # try:
    #     for k, v in data["parts"].items():
    #         if 851383 in v:
    #             echo(k, v)
    #             echo(data["s_between"])
    #             echo("s within", data["s_within"])
    #             echo("here1")
    #             echo(len(data["parts"]))
    #             echo(data["parts"])
    #             quit()
    # except:
    #     pass
    n_parts = len(data["parts"])
    n_reads = len(data["reads"])
    if n_parts >= 2:
        # Processed single edges and break apart connected
        for event in multi(data, bam, insert_size, insert_stdev, clip_length, min_support):
            yield event

    elif n_parts == 1:
        # Single isolated node, multiple reads

        d = data["parts"][next(iter(data["parts"]))]  # Get first key

        # rds = {}
        to_collect = {}
        for v in d:
            if v not in data["reads"]:
            # try:
            #     rds[v] = data["reads"][v]  # May have been collected already

            # except KeyError:
                to_collect[v] = data["n2n"][v]

        # rds.update(cy_graph.get_reads(bam, to_collect))
        # data["reads"] = rds
        data["reads"].update(graph.get_reads(bam, to_collect))
        yield single(bam, data, insert_size, insert_stdev, clip_length, min_support)
    #
    #
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
        else:
            reads_10kb = 0
    else:
        # Calculate max
        chrom_i = infile.get_tid(r["chrA"])
        if chrom_i in regions_depth:
            reads_10kb_left = coverage.calculate_coverage(r["posA"] - 10000, r["posA"] + 10000, regions_depth[chrom_i])
        else:
            reads_10kb_left = 0

        chrom_i = infile.get_tid(r["chrB"])
        if chrom_i in regions_depth:
            reads_10kb_right = coverage.calculate_coverage(r["posB"] - 10000, r["posB"] + 10000, regions_depth[chrom_i])
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
