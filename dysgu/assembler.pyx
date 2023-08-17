#distutils: language = c++
#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True
"""
A basic assembler/consensus sequence generator. Takes an overlap graph and merges reads in-place in a POA style.
"""

import warnings
import array

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

from dysgu.map_set_utils import timeit, echo

from dysgu.scikitbio._ssw_wrapper import StripedSmithWaterman
from libcpp.vector cimport vector as cpp_vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair

from libc.math cimport exp
from libc.stdlib cimport abs as c_abs
from libc.stdint cimport uint32_t, int32_t, uint64_t

from dysgu cimport map_set_utils
from dysgu.map_set_utils cimport DiGraph, unordered_set, unordered_map, EventResult
from dysgu.map_set_utils cimport hash as xxhasher

from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_seqi, bam_get_seq, bam_get_cigar

import numpy as np

ctypedef cpp_vector[int] int_vec_t

ctypedef map_set_utils.Py_Int2IntMap Py_Int2IntMap
ctypedef map_set_utils.Py_IntSet Py_IntSet

ctypedef EventResult EventResult_t


cdef extern from "graph_objects.hpp" nogil:

    cdef cppclass TwoWayMap:
        TwoWayMap() nogil
        uint64_t key_2_64(char, uint64_t, uint64_t, uint64_t) nogil
        void insert_tuple_key(uint64_t, int) nogil
        int has_tuple_key(uint64_t) nogil
        int get_index_prev() nogil
        uint64_t get_key_prev() nogil
        void idx_2_vec(int, cpp_vector[int]) nogil
        void key_2_vec(uint64_t, cpp_vector[int]) nogil

    cdef void graph_node_2_vec(uint64_t, cpp_vector[int]) nogil


basemap = np.array([ '.', 'A', 'C', '.', 'G', '.', '.', '.', 'T', '.', '.', '.', '.', '.', 'N', 'N', 'N'])
lowermap = np.array([ '.', 'a', 'c', '.', 'g', '.', '.', '.', 't', '.', '.', '.', '.', '.', 'n', 'n', 'n'])


cdef trim_cigar(uint32_t cigar_l, uint32_t *cigar_p, int pos, int approx_pos):
    cdef int seq_index = 0
    cdef int seq_start = 0
    cdef int index = 0
    cdef int start_index = 0
    cdef uint32_t cigar_value
    cdef int start_pos = pos
    cdef int end_index = cigar_l - 1
    cdef bint started = False
    cdef int opp, length
    aligned, large_gaps, n_small_gaps = 0, 0, 0
    for i in range(cigar_l):
        cigar_value = cigar_p[i]
        opp = <int> cigar_value & 15
        length = <int> cigar_value >> 4
        if opp == 4 and index == 0:
            seq_index += length
        if opp == 1:
            seq_index += length
        elif opp == 2:
            if not started:
                start_pos += length
            pos += length
        elif opp == 0 or opp == 7 or opp == 8 or opp == 3:
            if not started:
                if abs(pos + length - approx_pos) < 500:
                    started = True
                    start_pos = pos
                    pos += length
                    start_index = index
                    seq_start = seq_index
                else:
                    pos += length
                seq_index += length
            else:
                if abs(pos + length - approx_pos) > 500:
                    end_index = index + 1
                    break
                else:
                    pos += length
                    seq_index += length
        index += 1
    return start_index, end_index, start_pos, seq_start


cdef void add_to_graph(DiGraph& G, AlignedSegment r, cpp_vector[int]& nweight, TwoWayMap& ndict_r2,
                       int32_t approx_position, int32_t max_distance):

    cdef int i = 0

    cdef char* char_ptr_rseq = <char*>bam_get_seq(r._delegate)  # without cast, compiler complains of unsigned char*
    cdef char base

    cdef const unsigned char[:] quals = r.query_qualities

    cdef list cigar = r.cigartuples
    cdef int pos = r.pos
    cdef uint64_t current_pos = pos + 1
    cdef uint64_t o, key

    cdef int begin, end, step, p, qual, opp, length, n
    cdef bint start = True

    cdef int prev_node = -1
    cdef int ref_bases = 0

    cdef int target_bases = 2*max_distance
    cdef str seq
    cdef tuple k
    cdef bint done = 0
    cdef cpp_vector[int] vv = [0, 0, 0, 0]

    cdef int r_end = r.reference_end
    if r.reference_end + 100 < approx_position or approx_position < current_pos - 100:
        return  # shouldn't happen

    cdef uint32_t cigar_value
    cdef uint32_t cigar_l = r._delegate.core.n_cigar
    cdef uint32_t *cigar_p = bam_get_cigar(r._delegate)
    cdef int cigar_start = 0
    cdef int cigar_end = cigar_l
    if approx_position - pos > 500:
        cigar_start, cigar_end, current_pos, i = trim_cigar(cigar_l, cigar_p, pos, approx_position)
    cdef int cigar_index

    for cigar_index in range(cigar_start, cigar_end):
        cigar_value = cigar_p[cigar_index]
        opp = <int> cigar_value & 15
        length = <int> cigar_value >> 4
        if done:
            break
        if opp == 4:
            if start:
                for o in range(length, 0, -1):
                    qual = quals[i]
                    base = bam_seqi(char_ptr_rseq, i)
                    i += 1
                    # 0 = left soft clip
                    key = ndict_r2.key_2_64(base, current_pos, o, <uint64_t>0)
                    if ndict_r2.has_tuple_key(key):
                        n = ndict_r2.get_index_prev()
                    else:
                        n = G.addNode()
                        if n >= nweight.size():
                            nweight.push_back(0)
                        ndict_r2.insert_tuple_key(key, n)
                    nweight[n] += qual
                    if prev_node != -1:
                        G.updateEdge(prev_node, n, qual)
                    prev_node = n

            else:
                for o in range(1, length + 1, 1):
                    qual = quals[i]
                    base = bam_seqi(char_ptr_rseq, i)
                    i += 1
                    # 1 = right soft clip
                    key = ndict_r2.key_2_64(base, current_pos, o, <uint64_t>1)
                    if ndict_r2.has_tuple_key(key):
                        n = ndict_r2.get_index_prev()
                    else:
                        n = G.addNode()
                        if n >= nweight.size():
                            nweight.push_back(0)
                        ndict_r2.insert_tuple_key(key, n)
                    nweight[n] += qual
                    if prev_node != -1:
                        G.updateEdge(prev_node, n, qual)
                    prev_node = n

        elif opp == 1:  # Insertion
            if c_abs(<int32_t>current_pos - approx_position) > max_distance:
                i += length
                if current_pos > approx_position:
                    break  # out of range
                continue

            for o in range(1, length + 1, 1):
                qual = quals[i]
                base = bam_seqi(char_ptr_rseq, i)
                i += 1
                # 2 = insertion
                key = ndict_r2.key_2_64(base, current_pos, o, <uint64_t>2)
                if ndict_r2.has_tuple_key(key):
                    n = ndict_r2.get_index_prev()
                else:
                    n = G.addNode()
                    if n >= nweight.size():
                        nweight.push_back(0)
                    ndict_r2.insert_tuple_key(key, n)
                nweight[n] += qual
                if prev_node != -1:
                    G.updateEdge(prev_node, n, qual)
                prev_node = n
            # current_pos += 1  # <-- Reference pos increases 1

        elif opp == 2: # deletion
            current_pos += length # + 1

        elif opp == 0 or opp == 7 or opp == 8 or opp == 3:  # All match, match (=), mis-match (X), N's
            if current_pos < approx_position and current_pos + length < approx_position - max_distance: # abs(<int32_t>current_pos - approx_position + length) > max_distance:
                i += length
                current_pos += length
                continue

            for p in range(current_pos, current_pos + length):
                current_pos = p
                if current_pos < approx_position and approx_position - current_pos > max_distance:
                    i += 1
                    continue
                elif current_pos > approx_position and current_pos - approx_position > max_distance:
                    break
                ref_bases += 1
                if ref_bases > target_bases:
                    done = 1
                    break
                qual = quals[i]
                base = bam_seqi(char_ptr_rseq, i)
                i += 1
                key = ndict_r2.key_2_64(base, current_pos, <uint64_t>0, <uint64_t>4)
                if ndict_r2.has_tuple_key(key):
                    n = ndict_r2.get_index_prev()
                else:
                    n = G.addNode()
                    if n >= nweight.size():
                        nweight.push_back(0)
                    ndict_r2.insert_tuple_key(key, n)
                nweight[n] += qual
                if prev_node != -1:
                    G.updateEdge(prev_node, n, qual)
                prev_node = n

            current_pos += 1

        start = False


cdef int topo_sort2(DiGraph& G, cpp_deque[int]& order): #  except -1:

    cdef unordered_set[int] seen
    cdef unordered_set[int] explored

    # cdef cpp_deque[int] order
    cdef cpp_vector[int] fringe
    cdef cpp_vector[int] new_nodes
    cdef cpp_vector[int] neighbors
    cdef int v, n, w

    cdef cpp_vector[int] debug_res

    # with nogil:

    for v in range(G.numberOfNodes()):  # process all vertices in G
        if explored.find(v) != explored.end():
            continue

        fringe.clear()
        fringe.push_back(v)  # nodes yet to look at

        while fringe.size() != 0:

            w = fringe.back() # depth first search
            if explored.find(w) != explored.end():  # already looked down this branch
                fringe.pop_back()
                continue

            seen.insert(w)

            # Check successors for cycles and for new nodes
            if new_nodes.size() > 0:
                new_nodes.clear()

            G.neighbors(w, neighbors)
            for n in neighbors:
                if explored.find(n) == explored.end():

                    if seen.find(n) != seen.end(): #CYCLE !!
                        order.clear()
                        order.push_back(-1)
                        order.push_back(n)
                        order.push_back(w)
                        # return order
                        graph_node_2_vec(n, debug_res)
                        raise ValueError("Graph contains a cycle. Please report this. n={}, w={}, v={}. Node info n was: {}, {}, {}, {}".format(n, w, v, debug_res[0], debug_res[1], debug_res[2], debug_res[4]))

                    new_nodes.push_back(n)

            if new_nodes.size() > 0:  # Add new_nodes to fringe
                fringe.insert(fringe.end(), new_nodes.begin(), new_nodes.end())  # Extend

            else:  # No new nodes so w is fully explored
                explored.insert(w)

                order.push_front(w)
                fringe.pop_back()  # done considering this node


cdef cpp_deque[int] score_best_path(DiGraph& G, cpp_deque[int]& nodes_to_visit, cpp_vector[int]& n_weights):

    cdef cpp_vector[int] node_scores = n_weights  # Copy; total base quality at each node
    cdef Py_Int2IntMap pred_trace2 = map_set_utils.Py_Int2IntMap()

    cdef int best_score = -1
    cdef int best_node = -1
    cdef int i, u, node_weight, maxi, score

    cdef cpp_vector[cpp_pair[int, int]] neighborList  # node, weight
    cdef int pred_score, local_score, best_local_score, best_local_i
    cdef cpp_pair[int, int] pred
    cdef cpp_deque[int] path

    cdef int len_nodes = nodes_to_visit.size()

    if len_nodes == 0:
        return path

    with nogil:
        for i in range(0, len_nodes):
            u = nodes_to_visit[i]
            node_weight = n_weights[u]

            # Find best incoming node scores, best inEdge, and also best predecessor node
            G.forInEdgesOf(u, neighborList)
            if neighborList.size() == 0:
                node_scores[u] = node_weight
                if node_weight >= best_score:
                    best_score = node_weight
                    best_node = u
                continue

            best_local_score = -1
            best_local_i = -1
            for pred in neighborList:

                pred_score = node_scores[pred.first]  # first is other end of edge, second is weight
                score = node_weight + pred_score
                if score > node_scores[u]:
                    node_scores[u] = score

                if score >= best_score:
                    best_score = score
                    best_node = u

                # The sum of base-qualities of previous node
                local_score = pred.second
                if local_score > best_local_score:
                    best_local_score = local_score
                    best_local_i = pred.first

            if best_local_i != -1:
                pred_trace2.insert(u, best_local_i)

    if best_node == -1:
        return path
    # Start traceback from best scoring node, use locally best edge to trace path back
    u = best_node
    while True:
        path.push_front(u)
        if pred_trace2.has_key(u) == 0:
            break
        u = pred_trace2.get(u)

    return path


cdef dict get_consensus(rd, int position, int max_distance):

    cdef str seq = ""
    cdef int ref_start = -1
    cdef int ref_end = -1
    cdef int longest_left_sc, longest_right_sc
    cdef int begin = 0
    cdef int return_code

    cdef DiGraph G = DiGraph()
    cdef TwoWayMap ndict_r2
    cdef cpp_vector[int] node_weights

    for r in rd:
        if r.seq is None:
            continue
        if r.query_qualities is None or len(r.seq) != len(r.query_qualities):
            r.query_qualities = array.array("B", [1] * len(r.seq))

        add_to_graph(G, r, node_weights, ndict_r2, position, max_distance)

    cdef cpp_deque[int] nodes_to_visit2

    return_code = topo_sort2(G, nodes_to_visit2)

    if return_code == -1 or nodes_to_visit2.size() < 50:
        return {}

    cdef cpp_deque[int] path2

    path2 = score_best_path(G, nodes_to_visit2, node_weights)

    if path2.size() < 50:
        return {}

    cdef int front = path2.front()
    cdef int back = path2.back()

    # vec has the form base, current_pos, offset, base-type (left-clip/right-clip/insertion)
    cdef cpp_vector[int] vec = [0, 0, 0, 0]
    ndict_r2.idx_2_vec(front, vec)

    longest_left_sc = vec[2]
    vec.assign(vec.size(), 0)

    ndict_r2.idx_2_vec(back, vec)

    longest_right_sc = vec[2]
    vec.assign(vec.size(), 0)

    cdef tuple t
    cdef int item

    cdef str sequence = ""
    cdef int m, u, w

    cdef int count = 0
    cdef int finish = path2.size()

    cdef cpp_vector[float] path_qual = [1] * path2.size()

    for item in path2:

        ndict_r2.idx_2_vec(item, vec)

        if count == 0:
            u = -1
        else:
            u = path2[count - 1]
        if count == finish - 1:
            w = -1
        else:
            w = path2[count + 1]

        path_qual[count] = G.node_path_quality(u, item, w)

        if vec[3] != 4:
            m = vec[0]
            sequence += lowermap[m]
        else:
            m = vec[0]
            if ref_start == -1:
                ref_start = vec[1]

            elif vec[1] > ref_end:
                ref_end = vec[1]

            sequence += basemap[m]

        vec.assign(vec.size(), 0)

        count += 1

    seq = sequence

    # Trim off bad sequence
    cdef int i, start_seq, end_seq
    cdef int original_right_sc = longest_right_sc
    end_seq = len(seq)

    if longest_right_sc > 0:
        for i in range(len(seq) - 1, len(seq) - longest_right_sc, -1):
            if path_qual[i] < 0.5:
                end_seq = i
        end_seq = min(len(seq) - longest_right_sc + 300, end_seq)
        longest_right_sc -= len(seq) - end_seq

    start_seq = 0
    if longest_left_sc > 0:
        for i in range(longest_left_sc):
            if path_qual[i] < 0.5:
                start_seq = i
        start_seq = max(start_seq, longest_left_sc - 300)
        longest_left_sc -= start_seq

    # Average quality weight over soft-clip portion
    cdef float left_clip_weight = 0
    if longest_left_sc > 0:
        for i in range(start_seq, start_seq + longest_left_sc):
            left_clip_weight += node_weights[i]
        left_clip_weight = left_clip_weight / longest_left_sc

    cdef float right_clip_weight = 0
    if longest_right_sc > 0:
        for i in range(len(seq) - original_right_sc, end_seq):
            right_clip_weight += node_weights[i]
        right_clip_weight = right_clip_weight / longest_right_sc

    if start_seq != 0 or end_seq != len(seq):
        seq = seq[start_seq:end_seq]

    return {"contig": seq,
            "left_clips": longest_left_sc,
            "right_clips": longest_right_sc,
            "ref_bases": len(seq) - longest_left_sc - longest_right_sc,
            "ref_start": ref_start,
            "ref_end": ref_end,
            "bamrname": rd[0].rname,
            "left_weight": left_clip_weight,
            "right_weight": right_clip_weight}


cdef trim_sequence_from_cigar(r, int approx_pos, int max_distance):

    ct = r.cigartuples
    seq = r.seq
    cdef int seq_index = 0
    cdef int seq_start = 0

    cdef int index = 0  # current cigar index
    cdef int start_index = 0  # index into cigar

    cdef int start_pos = r.pos  # starting genome position of cigar at start_index
    cdef int pos = start_pos  # current genome position
    cdef int end_index = len(ct) - 1
    cdef bint started = False
    cdef int opp, length

    cdef int longest_left_sc = 0
    cdef int longest_right_sc = 0

    parts = []
    for opp, length in ct:

        if opp == 4 and index == 0:
            if abs(pos - approx_pos) < 50:
                # include clipped but trim
                longest_left_sc = 150 if length >= 150 else length
                parts.append(seq[length - 150 if length - 150 > 0 else 0: length])
            seq_index += length
            index += 1
            continue
        if opp == 4 and index > 0:
            longest_right_sc = 150 if length >= 150 else length
            e = 150 if length > 150 else length
            parts.append(seq[seq_index: seq_index + e])
            break

        if opp == 1:
            if started:
                parts.append(r.seq[seq_index:seq_index + length].lower())
            seq_index += length

        elif opp == 2:
            if not started:
                start_pos += length
            pos += length

        elif opp == 0 or opp == 7 or opp == 8 or opp == 3:

            if not started:
                if abs(pos + length - approx_pos) < 300:
                    parts.append(r.seq[seq_index:seq_index + length])  # upper implied

                    started = True
                    start_pos = pos
                    pos += length
                    start_index = index
                    seq_start = seq_index

                else:
                    pos += length

                seq_index += length

            else:
                if abs(pos + length - approx_pos) > 300:
                    end_index = index + 1
                    break
                else:
                    parts.append(r.seq[seq_index:seq_index + length])
                    pos += length
                    seq_index += length

        index += 1

    return {"contig": "".join(parts),
            "left_clips": longest_left_sc,
            "right_clips": longest_right_sc,
            "ref_bases": len(seq) - longest_left_sc - longest_right_sc,
            "ref_start": r.pos,
            "ref_end": r.reference_end,
            "bamrname": r.rname,
            "left_weight": 0,
            "right_weight": 0}


cpdef dict base_assemble(rd, int position, int max_distance):
    if len(rd) == 1:
        r = rd[0]
        rseq = r.seq
        ct = r.cigartuples
        if rseq is None or ct is None:
            return {}
        if len(rseq) > 250:
            return trim_sequence_from_cigar(r, position, max_distance)
        else:
            longest_left_sc = 0
            longest_right_sc = 0
            seq = ""
            begin = 0
            for opp, length in ct:
                if opp == 4 or opp == 1:
                    seq += rseq[begin:begin + length].lower()
                    if opp == 4:
                        if begin == 0:
                            longest_left_sc = length
                        else:
                            longest_right_sc = length
                    begin += length
                elif opp == 0 or opp == 7 or opp == 8 or opp == 3:
                    seq += rseq[begin:begin + length]
                    begin += length
            return {"contig": seq,
                    "left_clips": longest_left_sc,
                    "right_clips": longest_right_sc,
                    "ref_bases": len(seq) - longest_left_sc - longest_right_sc,
                    "ref_start": r.pos,
                    "ref_end": r.reference_end,
                    "bamrname": r.rname,
                    "left_weight": 0,
                    "right_weight": 0}
    return get_consensus(rd, position, max_distance)


cpdef float compute_rep(seq):

    cdef unordered_map[float, int] last_visited
    cdef float tot_amount = 0
    cdef float total_seen = 0
    cdef int k, i, diff
    cdef float decay, max_amount, amount

    cdef bytes s_bytes = bytes(seq.encode("ascii"))
    cdef const unsigned char* sub_ptr = s_bytes

    for k in (2, 3, 4, 5, 6):

        decay = 0.25 * 1/k
        max_amount = exp(-decay) * k  # If last kmer was the same as current kmer

        sub_ptr = s_bytes
        for i in range(len(seq) - k):

            a = xxhasher(sub_ptr, k, 42)
            if last_visited.find(a) != last_visited.end():
                diff = i - last_visited[a]
                x = exp(-decay * diff)
                amount = (k * x) / max_amount

            else:
                amount = 0
            if i > k:
                tot_amount += amount
                total_seen += 1
            last_visited[a] = i
            sub_ptr += 1

    if total_seen == 0:
        return 0

    return tot_amount / total_seen

cdef tuple get_rep(contig_seq):

    cdef int left_clip_end = 0
    cdef int right_clip_start = 0
    cdef float aligned_portion, clip_rep, clip_seen
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

    if clip_seen != 0:
        clip_rep = clip_rep / clip_seen
    return aligned_portion, clip_rep, right_clip_start - left_clip_end


def contig_info(events):

    cdef EventResult_t e
    for i in range(len(events)):
        e = events[i]
        gc_count = 0
        seq_length = 0
        if e.contig:
            cont = e.contig.upper()
            seq_length += len(cont)
            for letter in cont:
                if letter == "G" or letter == "C":
                    gc_count += 1

        if e.contig2:
            cont = e.contig2.upper()
            seq_length += len(cont)
            for letter in cont:
                if letter == "G" or letter == "C":
                    gc_count += 1

        if seq_length > 0:
            e.gc = round((gc_count / seq_length) * 100, 2)
        else:
            e.gc = 0

        sc_rep = 0
        aln_rep = 0
        aligned = 0
        seen = 0
        if e.contig:
            aln_rep1, sc_rep1, aligned_bases = get_rep(e.contig)
            sc_rep += sc_rep1
            aln_rep += aln_rep1
            aligned += aligned_bases
            seen += 1

        if e.contig2:
            aln_rep1, sc_rep1, aligned_bases = get_rep(e.contig2)
            seen += 1
            sc_rep += sc_rep1
            aln_rep += aln_rep1
            aligned += aligned_bases

        if seen > 0:
            aln_rep = aln_rep / seen
            sc_rep = sc_rep / seen

        e.rep = round(aln_rep, 3)
        e.rep_sc = round(sc_rep, 3)
        e.ref_bases = aligned

    return events


def check_contig_match(a, b, rel_diffs=False, diffs=8, ol_length=21, supress_seq=True, return_int=False):
    if not a or not b or len(a) > 5000 or len(b) > 5000:
        return 0
    query = StripedSmithWaterman(str(a), suppress_sequences=supress_seq)
    alignment = query(str(b))

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

    rel_diff = (aln_s - total_overhangs) / expected
    diff = expected - aln_s + total_overhangs

    if aln_s > 100 and total_overhangs < qe - qs:
        if rel_diff > 0.7:
            if return_int:
                return 1
            return (qs, qe, als, ale, alignment.cigar, alignment.aligned_query_sequence,
                    alignment.aligned_target_sequence)
        return 0

    if diff > diffs:  # e.g. 2 mis-matches + 2 unaligned overhanging bits
        return 0
    else:
        if return_int:
            return 1
        return (qs, qe, als, ale,
                alignment.cigar,
                alignment.aligned_query_sequence,
                alignment.aligned_target_sequence)

