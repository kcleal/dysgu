#distutils: language = c++
#cython: language_level = 3
"""
A basic assembler. Takes an overlap graph and merges reads in-place in a POA style. Different soft-clipped regions
are then overlapped and 'linked'.
"""
import networkit as nk

import numpy as np
cimport numpy as np
from collections import deque
import click
from skbio.alignment import StripedSmithWaterman
from libcpp.vector cimport vector as cpp_vector
# from preshed.maps cimport PreshMap

from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair

from libc.math cimport exp


def echo(*args):
    click.echo(args, err=True)

DTYPE = np.int64
ctypedef np.int64_t DTYPE_t
ctypedef cpp_vector[int] int_vec_t
ctypedef cpp_pair[int, int] get_val_result
ctypedef cpp_pair[int, int] cpp_item

# ctypedef Py_Int2IntVecMap[int, int_vec_t] node_dict_t
# ctypedef Py_IntVec2IntMap[int_vec_t, int] node_dict2_r_t



cdef extern from "wrap_map_set2.h":
    cdef cppclass Int2IntMap:
        Int2IntMap()
        void insert(int, int)
        void erase(int)
        int has_key(int)
        int get(int)
        get_val_result get_value(int key)
        int size()


cdef class Py_Int2IntMap:
    """Fast 32bit integer to 32bit integer unordered map using tsl::robin-map"""
    cdef Int2IntMap *thisptr
    def __cinit__(self):
        self.thisptr = new Int2IntMap()
    def __dealloc__(self):
        del self.thisptr
    cpdef void insert(self, int key, int value):
        self.thisptr.insert(key, value)
    cpdef void erase(self, int key):
        self.thisptr.erase(key)
    cpdef int has_key(self, int key):
        return self.thisptr.has_key(key)
    cpdef int get(self, int key):
        return self.thisptr.get(key)
    cpdef get_val_result get_value(self, int key):
        return self.thisptr.get_value(key)
    cpdef int size(self):
        return self.thisptr.size()


cdef extern from "wrap_map_set2.h":
    cdef cppclass IntSet:
        IntSet()
        void insert(int)
        void erase(int)
        int has_key(int)
        int get(int)
        int size()


cdef class Py_IntSet:
    """Fast 32 bit int set using tsl::robin-set"""
    cdef IntSet *thisptr
    def __cinit__(self):
        self.thisptr = new IntSet()
    def __dealloc__(self):
        del self.thisptr
    cpdef void insert(self, int key):
        self.thisptr.insert(key)
    cpdef void erase(self, int key):
        self.thisptr.erase(key)
    cpdef int has_key(self, int key):
        return self.thisptr.has_key(key)
    cpdef int size(self):
        return self.thisptr.size()

# Vector as key doesnt work so well, hashing is too slow
# cdef extern from "wrap_map_set2.h":
#     cdef cppclass Int2IntVecMap:
#         Int2IntVecMap()
#         void insert(int, cpp_vector[int])
#         void erase(int)
#         int has_key(int)
#         cpp_vector[int] get(int)
#         int size()
#
# cdef class Py_Int2IntVecMap:
#     """Fast 32 bit int key to 32 bit vector value using tsl::robin-map"""
#     cdef Int2IntVecMap *thisptr
#     def __cinit__(self):
#         self.thisptr = new Int2IntVecMap()
#     def __dealloc__(self):
#         del self.thisptr
#     cpdef void insert(self, int key, cpp_vector[int] value):
#         self.thisptr.insert(key, value)
#     cpdef void erase(self, int key):
#         self.thisptr.erase(key)
#     cpdef int has_key(self, int key):
#         return self.thisptr.has_key(key)
#     cpdef cpp_vector[int] get(self, int key):
#         return self.thisptr.get(key)
#     cpdef int size(self):
#         return self.thisptr.size()
#
#
# cdef extern from "wrap_map_set2.h":
#     cdef cppclass IntVec2IntMap:
#         IntVec2IntMap()
#         void insert(cpp_vector[int], int)
#         void erase(cpp_vector[int])
#         int has_key(cpp_vector[int])
#         int get(cpp_vector[int])
#         int size()
#
# cdef class Py_IntVec2IntMap:
#     """Fast 32 bit int vec key to 32 bit value using tsl::robin-map"""
#     cdef IntVec2IntMap *thisptr
#     def __cinit__(self):
#         self.thisptr = new IntVec2IntMap()
#     def __dealloc__(self):
#         del self.thisptr
#     cpdef void insert(self, cpp_vector[int] key, int value):
#         self.thisptr.insert(key, value)
#     cpdef void erase(self, cpp_vector[int] key):
#         self.thisptr.erase(key)
#     cpdef int has_key(self, cpp_vector[int] key):
#         return self.thisptr.has_key(key)
#     cpdef int get(self, cpp_vector[int] key):
#         return self.thisptr.get(key)
#     cpdef int size(self):
#         return self.thisptr.size()


cpdef void add_to_graph(G, r, cpp_vector[int]& nweight, ndict_r):

    cdef int i = 0
    cdef str rseq = r.seq
    cdef const unsigned char[:] quals = r.query_qualities

    cdef list cigar = r.cigartuples
    cdef int pos = r.pos
    cdef int current_pos = pos + 1

    cdef int o, begin, end, step, p, qual, opp, length
    cdef int start = 1

    cdef int prev_node = -1

    cdef str seq  # left clip == 0, right clip == 1, insertion == 2, match == 4
    cdef tuple k

    for opp, length in cigar:

        if opp == 4:
            if start:
                for o in range(length, 0, -1):
                    seq = rseq[i]
                    qual = quals[i]
                    i += 1

                    k = (seq, current_pos, o, 0)

                    if k in ndict_r:  #
                        n = ndict_r[k]  #

                    else:
                        n = G.addNode()
                        if n >= nweight.size():
                            nweight.push_back(0)

                        # ndict[n] = k  #
                        ndict_r[k] = n  #

                    nweight[n] += qual
                    if prev_node != -1:
                        G.addEdge(prev_node, n)

                    prev_node = n

            else:
                for o in range(1, length + 1, 1):
                    seq = rseq[i]
                    qual = quals[i]
                    i += 1

                    k = (seq, current_pos, o, 1)

                    if k in ndict_r:
                        n = ndict_r[k]

                    else:
                        n = G.addNode()
                        if n >= nweight.size():
                            nweight.push_back(0)

                        # ndict[n] = k
                        ndict_r[k] = n

                    nweight[n] += qual
                    if prev_node != -1:
                        G.addEdge(prev_node, n)

                    prev_node = n

        elif opp == 1:  # Insertion

            for o in range(1, length + 1, 1):

                seq = rseq[i]
                qual = quals[i]
                i += 1

                k = (seq, current_pos, o, 2)

                if k in ndict_r:
                    n = ndict_r[k]

                else:
                    n = G.addNode()
                    if n >= nweight.size():
                        nweight.push_back(0)

                    # ndict[n] = k
                    ndict_r[k] = n

                nweight[n] += qual
                if prev_node != -1:
                    G.addEdge(prev_node, n)
                prev_node = n

            current_pos += 1  # Reference pos increases only 1

        elif opp == 2 or opp == 5:  # Hard clip or deletion
            current_pos += length

        elif opp == 0 or opp == 7 or opp == 8 or opp == 3:  # All match, match (=), mis-match (X), N's

            for p in range(current_pos, current_pos + length):

                current_pos = p
                seq = rseq[i]
                qual = quals[i]
                i += 1

                k = (seq, current_pos, 0, 4)

                if k in ndict_r:
                    n = ndict_r[k]

                else:
                    n = G.addNode()
                    if n >= nweight.size():
                        nweight.push_back(0)

                    # ndict[n] = k
                    ndict_r[k] = n
                    #

                nweight[n] += qual
                if prev_node != -1:
                    G.addEdge(prev_node, n)

                prev_node = n

        start = 0

cdef cpp_deque[int] topo_sort2(G):
    # https://networkx.github.io/documentation/networkx-1.9/_modules/networkx/algorithms/dag.html#topological_sort

    cdef Py_IntSet seen = Py_IntSet()
    cdef Py_IntSet explored = Py_IntSet()
    cdef cpp_deque[int] order
    cdef cpp_vector[int] fringe
    cdef cpp_vector[int] new_nodes
    cdef int v, n, w

    for v in G.nodes():     # process all vertices in G
        if explored.has_key(v) == 1:
            continue
        fringe.clear()
        fringe.push_back(v)   # nodes yet to look at
        while fringe.size() != 0:

            w = fringe.back()  # depth first search
            if explored.has_key(w) == 1: # already looked down this branch
                fringe.pop_back()
                continue
            seen.insert(w)     # mark as seen

            # Check successors for cycles and for new nodes
            if new_nodes.size() > 0:
                new_nodes.clear()

            for n in G.neighbors(w):
                if explored.has_key(n) == 0:
                    if seen.has_key(n) == 1: #CYCLE !!
                        raise ValueError("Graph contains a cycle.")
                    new_nodes.push_back(n)

            if new_nodes.size() > 0:   # Add new_nodes to fringe
                fringe.insert(fringe.end(), new_nodes.begin(), new_nodes.end())  # Extend

            else:           # No new nodes so w is fully explored
                explored.insert(w)

                order.push_front(w)
                fringe.pop_back()    # done considering this node

    return order


def score_best_path(G, nodes_to_visit, cpp_vector[int]& n_weights):

    # copy
    cdef cpp_vector[int] node_scores = n_weights
    cdef Py_Int2IntMap pred_trace2 = Py_Int2IntMap()

    cdef int best_score = -1
    cdef int best_node = -1
    cdef int i, u, node_weight, maxi, score

    for i in range(0, len(nodes_to_visit)):

        u = nodes_to_visit[i]
        node_weight = n_weights[u]

        # Find best incoming node scores
        # Todo find a way to cythonize this lambda - not possible to use with cdef
        neighborList = []
        G.forInEdgesOf(u, lambda unode, pred, edgeweight, edgeid: neighborList.append((node_scores[pred], pred,
                                                                                       n_weights[pred])))

        if len(neighborList) == 0:
            node_scores[u] = node_weight
            if node_weight >= best_score:
                best_score = node_weight
                best_node = u
            continue

        # Sum of best path in graph
        maxi = neighborList.index(max(neighborList))
        score = node_weight + neighborList[maxi][0]
        node_scores[u] = score
        if score >= best_score:
            best_score = score
            best_node = u

        # Also track best weighted local node
        maxi = neighborList.index(max(neighborList, key=lambda x: x[2]))
        pred_trace2.insert(u, neighborList[maxi][1])

    if best_node == -1:
        return []
    # Start traceback from best scoring node, use locally best edge to trace path back
    path = deque([])

    u = best_node
    while True:
        path.appendleft(u)
        # path2.push_back(u)
        if pred_trace2.has_key(u) == 0:
            break
        u = pred_trace2.get(u)

    return path


cpdef dict base_assemble(rd):

    # Note supplementary are included in assembly; helps link regions
    # Get reads of interest

    cdef str seq = ""
    cdef int ref_start = -1
    cdef int ref_end = -1
    cdef int longest_left_sc, longest_right_sc
    cdef int begin = 0

    G = nk.Graph(weighted=False, directed=True)

    node_dict_r = {}

    # cdef Py_Int2IntVecMap node_dict2 = Py_Int2IntVecMap()
    # cdef Py_IntVec2IntMap node_dict2_r = Py_IntVec2IntMap()

    cdef cpp_vector[int] node_weights

    if len(rd) == 1:
        r = rd[0]
        rseq = r.seq
        ct = r.cigartuples

        longest_left_sc = 0
        longest_right_sc = 0
        if ct[0][0] != 4 and ct[-1][0] != 4:
            return {}
        if ct[0][0] == 4:
            begin = ct[0][1]
            seq += rseq[:begin].lower()
            longest_left_sc = ct[0][1]

        if ct[-1][0] == 4:
            end = len(rseq) - ct[-1][-1]
            seq += rseq[begin:end]
            seq += rseq[end:].lower()
            longest_right_sc = ct[-1][1]
        else:
            seq += rseq[begin:len(rseq)]

        return {"contig": seq,
                "left_clips": longest_left_sc,
                "right_clips": longest_right_sc,
                "ref_bases": len(seq) - longest_left_sc - longest_right_sc,
                "ref_start": r.pos,
                "ref_end": r.reference_end,
                "bamrname": r.rname}

    for r in rd:
        try:
            r.seq
        except:
            continue
        if r.seq is None:
            continue
        if len(r.seq) != len(r.query_qualities):
            continue

        add_to_graph(G, r, node_weights, node_dict_r)

    cdef cpp_deque[int] nodes_to_visit2 = topo_sort2(G)

    node_list = []
    cdef int i
    for i in nodes_to_visit2:
        node_list.append(i)

    path2 = score_best_path(G, node_list, node_weights)

    node_dict = list(node_dict_r.keys())


    if len(path2) == 0:
        return {}
    longest_left_sc = node_dict[path2[0]][2]
    longest_right_sc = node_dict[path2[-1]][2]
    if longest_left_sc == 0 and longest_right_sc == 0:
            return {}  # No soft-clips, so not overlapping a break


    cdef tuple t
    cdef int item
    for item in path2:

        t = node_dict[item]

        if t[3] != 4:
            seq += t[0].lower()
        else:
            if ref_start == -1:
                ref_start = t[1]
            elif t[1] > ref_end:
                ref_end = t[1]
            seq += t[0]

    return {"contig": seq,
            "left_clips": longest_left_sc,
            "right_clips": longest_right_sc,
            "ref_bases": len(seq) - longest_left_sc - longest_right_sc,
            "ref_start": ref_start,
            "ref_end": ref_end,
            "bamrname": rd[0].rname}


cdef float compute_rep(str seq):

    last_visited = {}
    tot = []

    cdef int k, i, diff
    cdef float decay, max_amount, amount
    cdef str a

    for k in (2, 3, 4, 5, 6):

        decay = 0.25 * 1/k
        max_amount = exp(-decay) * k  # If last kmer was the same as current kmer

        for i in range(len(seq) - k):
            a = seq[i:i+k]

            if a in last_visited:
                diff = i - last_visited[a]
                amount = (((diff * exp(-decay * diff)) / diff) * k) / max_amount
            else:
                amount = 0

            if i > k:
                tot.append(amount)

            last_visited[a] = i

    if len(tot) == 0:
        return 0
    return sum(tot) / len(tot)


cdef tuple get_rep(str contig_seq):

    cdef int left_clip_end, right_clip_start
    cdef float aligned_portion, clip_rep
    cdef str clip
    # Get left clip
    for left_clip_end in range(len(contig_seq)):
        if contig_seq[left_clip_end].isupper():
            break


    for right_clip_start in range(len(contig_seq) - 1, -1, -1):
        if contig_seq[right_clip_start].isupper():
            right_clip_start += 1
            break

    aligned_portion = compute_rep(contig_seq[left_clip_end: right_clip_start])

    clip_rep = 0
    clip_seen = 0
    if left_clip_end > 0:
        clip = contig_seq[:left_clip_end]
        clip_rep += compute_rep(clip)
        clip_seen += 1

    elif right_clip_start < len(contig_seq):
        clip = contig_seq[right_clip_start:]
        clip_rep += compute_rep(clip)
        clip_seen += 1

    clip_rep = clip_rep / clip_seen
    return aligned_portion, clip_rep, right_clip_start - left_clip_end



cpdef list contig_info(list events):

    cdef int i, aligned, seen, aligned_bases
    cdef float gc_count, seq_length
    # cdef float c, n_conts
    cdef float sc_rep, aln_rep, sc_rep1, aln_rep1
    cdef str letter
    for i in range(len(events)):
        e = events[i]
        gc_count = 0
        seq_length = 0
        if e["contig"]:
            seq_length += <float>len(e["contig"])
            for letter in e["contig"]:
                if letter == "G" or letter == "C" or letter == "c" or letter == "g":
                    gc_count += 1
        if e["contig2"]:
            seq_length += <float>len(e["contig2"])
            for letter in e["contig2"]:
                if letter == "G" or letter == "C" or letter == "c" or letter == "g":
                    gc_count += 1
        if seq_length > 0:
            e["gc"] = round((gc_count / seq_length) * 100, 2)
        else:
            e["gc"] = 0

        sc_rep = 0
        aln_rep = 0
        aligned = 0
        seen = 0
        if e["contig"]:
            aln_rep1, sc_rep1, aligned_bases = get_rep(e["contig"])
            sc_rep += sc_rep1
            aln_rep += aln_rep1
            aligned += aligned_bases
            seen += 1

        if e["contig2"]:
            aln_rep1, sc_rep1, aligned_bases = get_rep(e["contig2"])
            seen += 1
            sc_rep += sc_rep1
            aln_rep += aln_rep1
            aligned += aligned_bases

        aln_rep = aln_rep / seen
        sc_rep = sc_rep / seen

        e["rep"] = round(aln_rep, 3)
        e["rep_sc"] = round(sc_rep, 3)
        e["ref_bases"] = aligned

    return events


def check_contig_match(a, b, diffs=8, ol_length=21, supress_seq=True, return_int=False):

    query = StripedSmithWaterman(str(a), suppress_sequences=supress_seq)
    alignment = query(str(b))
    # echo(alignment)
    qs, qe = alignment.query_begin, alignment.query_end
    als, ale = alignment.target_begin, alignment.target_end_optimal

    # Find the length of any unaligned overhangs
    extent_left = min((qs, als))
    extent_right = min((len(a) - qe, len(b) - ale))
    total_overhangs = extent_left + extent_right
    aln_s = alignment.optimal_alignment_score
    expected = (qe - qs) * 2  # +2 is the score for a match

    if expected < 2 * ol_length:  # Match score * Minimum clip length
        return 0
    diff = expected - aln_s + total_overhangs  # old diff thresh = 8
    # diff = (expected - aln_s - total_overhangs) / expected
    if diff > diffs:  # e.g. 2 mis-matches + 2 unaligned overhanging bits
        return 0
    else:
        if return_int:
            return 1
        return (qs, qe, als, ale,
                alignment.cigar,
                alignment.aligned_query_sequence,
                alignment.aligned_target_sequence)


def get_upper_start_end(a):
    a_start, a_end = -1, 0
    for idx, l in enumerate(a):
        if l.isupper():
            if a_start == -1:
                a_start = idx
            if idx > a_end:
                a_end = idx
    return a_start, a_end + 1


def get_mark_result(res, insertion, a, b, sa, sb, a_start, a_end, b_start, b_end, b_rev=False):
    if insertion < 0:
        res["mark"] = "microh"
    else:
        res["mark"] = "ins"
    edit_dis = len([1 for i, j in zip(sa, sb) if i != j])
    res["mark_seq"] = sa
    res["mark_ed"] = edit_dis

    if insertion > 0:
        # Look for templated insertion
        a_align = a[a_start:a_end].replace("-", "")  # Remove any deletion markers
        b_align = b[b_start:b_end].replace("-", "")

        a_query = StripedSmithWaterman(sa)
        a_alignment = a_query(a_align)
        b_alignment = a_query(b_align)

        if a_alignment.optimal_alignment_score >= b_alignment.optimal_alignment_score:
            cont = "A"
            aln = a_alignment
        else:
            cont = "B"
            aln = b_alignment

        aqs = aln.aligned_query_sequence
        try:
            tqs = aln.aligned_target_sequence
        except IndexError:
            tqs = None  # No alignment

        if tqs:
            edit_dis = len([1 for i, j in zip(aqs, tqs) if i.upper() != j])
            if b_rev and cont == "B":

                v = "contB_rev:pos={}:ed={}:align={}".format(aln.target_begin, edit_dis, aqs)
            else:
                v = "cont{}:pos={}:ed={}:align={}".format(cont, aln.target_begin, edit_dis, aqs)
            l = aln.target_end_optimal - aln.target_begin + 1

            res["templated_ins_info"] = v
            res["templated_ins_len"] = l

    return res


def get_microh_or_ins(aln_idx):
    qs, qe, als, ale, q_cigar, q_aln, t_aln = aln_idx
    a = q_aln  # Use actual alignment sequence - keep deletions and insertions in place
    b = t_aln

    a_start, a_end = get_upper_start_end(a)
    b_start, b_end = get_upper_start_end(b)

    # Check for overlap of gap
    res = {"mark": "blunt", "mark_seq": "", "mark_ed": "", "templated_ins_info": "", "templated_ins_len": ""}
    if a_start >= b_start:
        insertion = a_start - b_end
        if insertion != 0:
            v = slice(*sorted([a_start, b_end]))
            sa = a[v]
            sb = b[v]
            res = get_mark_result(res, insertion, a, b, sa, sb, a_start, a_end, b_start, b_end)

    else:
        insertion = b_start - a_end
        if insertion != 0:
            v = slice(*sorted([b_start, a_end]))
            sa = a[v]
            sb = b[v]
            res = get_mark_result(res, insertion, a, b, sa, sb, a_start, a_end, b_start, b_end)

    return res


def link_pair_of_assemblies(a, b, clip_length):

    # Safest way is to try forward and reverse complement
    m = check_contig_match(a["contig"], b["contig"], supress_seq=False)

    if m != 0:
        a["linked"] = 1
        h = get_microh_or_ins(m)

    else:

        m = check_contig_match(a["contig"], b["contig_rev"], supress_seq=False)
        if m != 0:
            a["linked"] = 1
            h = get_microh_or_ins(m)
        else:
            a["linked"] = 0
            h = {"mark": "None", "mark_seq": "", "mark_ed": "", "templated_ins_info": "",
                 "templated_ins_len": ""}
    a.update(h)
    return a, b
