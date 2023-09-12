#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True
from __future__ import absolute_import
from collections import defaultdict, deque, namedtuple
import numpy as np
cimport numpy as np
import sortedcontainers
import cython
from cpython cimport array
import array
import re
import logging
from dysgu.map_set_utils cimport unordered_map as robin_map, Py_SimpleGraph
from dysgu.map_set_utils cimport multimap as cpp_map
from dysgu cimport map_set_utils
from dysgu.io_funcs import intersecter
from dysgu.map_set_utils cimport unordered_set, cigar_clip, clip_sizes_hard, is_reciprocal_overlapping, span_position_distance
from dysgu.map_set_utils cimport hash as xxhasher
from dysgu.map_set_utils cimport MinimizerTable
from dysgu.extra_metrics import BadClipCounter
from dysgu.map_set_utils import echo  # for debugging
from libcpp.string cimport string
from libcpp.deque cimport deque as cpp_deque
from libcpp.vector cimport vector
from libcpp.pair cimport pair as cpp_pair
from libcpp.unordered_map cimport unordered_map
from libc.stdint cimport uint16_t, uint32_t, int32_t, uint64_t
from libc.stdlib cimport abs as c_abs
from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_qname, bam_seqi, bam_get_seq, bam_get_cigar

ctypedef cpp_pair[int, int] cpp_item

ctypedef map_set_utils.Py_IntSet Py_IntSet
ctypedef map_set_utils.Py_Int2IntMap Py_Int2IntMap

ctypedef cpp_map[int, cpp_item] ScopeItem_t
ctypedef vector[ScopeItem_t] ForwardScope_t

ctypedef cpp_pair[int, cpp_item] event_item
ctypedef cpp_pair[long, int] cpp_long_pair
ctypedef long int long_int

ctypedef PairedEndScoper PairedEndScoper_t
ctypedef TemplateEdges TemplateEdges_t
ctypedef NodeToName NodeToName_t
ctypedef ClipScoper ClipScoper_t


# 1 = soft-clipped split-read, the supplementary mapping might be discarded. whole read is a node
# 2 = deletion event contained within read
# 3 = insertion event contained within read
# 4 = mate-pair is unmapped, or supplementary uninformative, probably an insertion
ctypedef enum ReadEnum_t:
    DISCORDANT = 0
    SPLIT = 1
    DELETION = 2
    INSERTION = 3
    BREAKEND = 4


cdef void sliding_window_minimum(int k, int m, str s, unordered_set[long]& found):
    """End minimizer
    https://github.com/keegancsmith/Sliding-Window-Minimum/blob/master/sliding_window_minimum.py"""
    cdef int i = 0
    cdef int end = len(s) - m + 1
    cdef int last_idx = end - 1
    cdef cpp_deque[cpp_long_pair] window2
    cdef long hx2
    cdef bytes s_bytes = bytes(s.encode("ascii"))
    cdef const unsigned char* sub_ptr = s_bytes
    with nogil:
        while i < end:
            hx2 = xxhasher(sub_ptr, m, 42)
            if i == 0 or i == last_idx:
                found.insert(hx2)
            while window2.size() != 0 and window2.back().first >= hx2:
                window2.pop_back()
            window2.push_back(cpp_long_pair(hx2, i))
            while window2.front().second <= i - k:
                window2.pop_front()
            sub_ptr += 1
            minimizer_i = window2.front().first
            found.insert(minimizer_i)
            i += 1


cdef str left_soft_clips(str seq, int code_length):
    return seq[0:code_length]


cdef str right_soft_clips(str seq, int code_length):
    return seq[len(seq) - code_length:]


cdef uint64_t pair_to_int64(uint32_t left, uint32_t right) nogil:
    return <uint64_t> left << 32 | right


cdef uint32_t int64_left(uint64_t item) nogil:
    cdef uint32_t v = 0
    v |= item >> 32
    return v

cdef uint32_t int64_right(uint64_t item) nogil:
    return <uint32_t> item


cdef class ClipScoper:
    # todo sort out types; use uint64_t instead of long
    # first in pair is chrom position, second is node name
    cdef cpp_deque[cpp_pair[int, int]] scope_left, scope_right
    # MinimizerTable - keys are read-names, value is a set of longs, minimizers associated with read
    cdef MinimizerTable read_minimizers_left, read_minimizers_right
    # Hashmap key is minimizer val, value is a set of chrom pos + node name packed as a 64 int
    cdef robin_map[uint64_t, unordered_set[uint64_t]] clip_table_left, clip_table_right
    cdef int max_dist, k, w, clip_length, current_chrom, minimizer_support_thresh, minimizer_matches
    cdef float target_density, upper_bound_n_minimizers
    cdef int n_local_minimizers_left, n_local_minimizers_right
    cdef object rds_clip
    cdef bint minimizer_mode

    """Keeps track of which reads are in scope"""
    def __init__(self, int max_dist, int k, int m, int clip_length, int minimizer_support_thresh,
                 int minimizer_breadth, int read_length):
        self.max_dist = max_dist
        self.k = k
        self.w = m
        self.clip_length = clip_length
        cdef long hx2
        cdef bytes s_bytes
        cdef const unsigned char* sub_ptr
        self.current_chrom = 0
        self.minimizer_support_thresh = minimizer_support_thresh
        self.minimizer_matches = minimizer_breadth
        self.n_local_minimizers_left = 0
        self.n_local_minimizers_right = 0
        self.target_density = 2. / (m + 1)
        self.upper_bound_n_minimizers = read_length * self.target_density
        self.rds_clip = [deque([]), deque([])]
        self.minimizer_mode = False

    cdef void _add_m_find_candidates(self, clip_seq, int name, int idx, int position,
                                     robin_map[uint64_t, unordered_set[uint64_t]]& clip_table,
                                     unordered_set[int]& clustered_nodes):
        cdef unordered_set[long] clip_minimizers
        cdef unordered_set[uint64_t].iterator targets_iter, targets_end
        sliding_window_minimum(self.k, self.w, clip_seq, clip_minimizers)
        if clip_minimizers.empty():
            return
        # add read minimizers from read, and find partners
        # idx 0 is for left clips, 1 is for right clips
        target_counts = defaultdict(int)
        cdef int total_m_found = 0
        cdef int find_candidate = 1
        cdef long m
        cdef long set_item
        cdef int n_local_minimizers, n_local_reads
        cdef float upper_bound
        if idx == 0:
            n_local_minimizers = self.n_local_minimizers_left
        else:
            n_local_minimizers = self.n_local_minimizers_right

        n_local_reads = self.scope_left.size() if idx == 0 else self.scope_right.size()
        upper_bound = (1 + (n_local_reads * 0.15)) * self.upper_bound_n_minimizers
        if n_local_minimizers > upper_bound:
            find_candidate = 0

        cdef cpp_pair[long, int] minimizer_pos_pair
        cdef int mitem, item_position
        for m in clip_minimizers:
            if idx == 0:
                self.read_minimizers_left.insert(name, m)
            else:
                self.read_minimizers_right.insert(name, m)
            if clip_table.find(m) == clip_table.end():
                clip_table[m].insert(pair_to_int64(position, name))
                if idx == 0:
                    self.n_local_minimizers_left += 1
                else:
                    self.n_local_minimizers_right += 1
                continue
            elif find_candidate:  # Look for suitable partners
                targets_iter = clip_table[m].begin()
                targets_end = clip_table[m].end()
                while targets_iter != targets_end:
                    set_item = dereference(targets_iter)
                    item_position = int64_left(set_item)
                    if abs(item_position - position) < 7:
                        mitem = int64_right(set_item)
                        total_m_found += 1
                        target_counts[mitem] += 1
                        support = (total_m_found / 2) + target_counts[mitem]
                        if support >= self.minimizer_support_thresh:
                            clustered_nodes.insert(mitem)
                        # Maximum edges for each read
                        if clustered_nodes.size() >= 5:
                            find_candidate = 0
                            break
                    preincrement(targets_iter)
            clip_table[m].insert(pair_to_int64(position, name))

    cdef void _refresh_scope(self, cpp_deque[cpp_pair[int, int]]& scope, int position, MinimizerTable& mm_table,
                             robin_map[uint64_t, unordered_set[uint64_t]]& clip_table, int index
                             ):
        # Remove out of scope reads and minimizers
        cdef int item_position, name
        cdef long m
        cdef cpp_pair[int, int] name_pair
        cdef unordered_set[long].iterator set_iter
        cdef unordered_set[long].iterator set_iter_end
        while True:
            if scope.empty():
                break
            if abs(scope[0].first - position) > self.max_dist:
                item_position = scope[0].first
                name = scope[0].second
                scope.pop_front()
                if mm_table.has_key(name):
                    set_iter = mm_table.get_iterator_begin()
                    set_iter_end = mm_table.get_iterator_end()
                    while set_iter != set_iter_end:
                        m = dereference(set_iter)
                        if clip_table.find(m) != clip_table.end():
                            clip_table[m].erase(pair_to_int64(item_position, name))
                            if clip_table[m].empty():
                                if index == 0:
                                    self.n_local_minimizers_left -= 1
                                else:
                                    self.n_local_minimizers_right -= 1
                        preincrement(set_iter)
                    mm_table.erase(name)
            else:
                break
    cdef void _insert(self, seq, int cigar_start, int cigar_end, int input_read, int position,
                      unordered_set[int]& clustered_nodes):
        # Find soft-clips of interest
        if cigar_start >= self.clip_length:
            clip_seq = left_soft_clips(seq, cigar_start)
            self._refresh_scope(self.scope_left, position, self.read_minimizers_left, self.clip_table_left, 0)
            self._add_m_find_candidates(clip_seq, input_read, 0, position, self.clip_table_left, clustered_nodes)
            self.scope_left.push_back(cpp_item(position, input_read))
        if cigar_end >= self.clip_length:
            clip_seq = right_soft_clips(seq, cigar_end)
            self._refresh_scope(self.scope_right, position, self.read_minimizers_right, self.clip_table_right, 1)
            self._add_m_find_candidates(clip_seq, input_read, 1, position, self.clip_table_right, clustered_nodes)
            self.scope_right.push_back(cpp_item(position, input_read))

    cdef void update(self, r, int input_read, int chrom, int position,
               unordered_set[int]& clustered_nodes):
        clip_left, clip_right = map_set_utils.clip_sizes(r)
        if chrom != self.current_chrom:
            self.scope_left.clear()
            self.scope_right.clear()
            self.clip_table_left.clear()
            self.clip_table_right.clear()
            self.n_local_minimizers_left = 0
            self.n_local_minimizers_right = 0
            self.current_chrom = chrom
        self._insert(r.seq, clip_left, clip_right, input_read, position, clustered_nodes)


cdef struct LocalVal:
    int chrom2
    int pos2
    int node_name
    int length_from_cigar
    ReadEnum_t read_enum


cdef LocalVal make_local_val(int chrom2, int pos2, int node_name, ReadEnum_t read_enum, int length_from_cigar) nogil:
    cdef LocalVal item
    item.chrom2 = chrom2
    item.pos2 = pos2
    item.node_name = node_name
    item.read_enum = read_enum
    item.length_from_cigar = length_from_cigar
    return item


cdef class PairedEndScoper:
    cdef int clst_dist
    cdef int max_dist
    cdef int local_chrom
    cdef cpp_map[int, LocalVal] loci  # Track the local breaks and mapping locations
    cdef vector[cpp_map[int, LocalVal]] chrom_scope  # Track the mate-pair breaks and locations
    cdef float norm
    cdef float thresh
    cdef bint paired_end
    def __init__(self, max_dist, clst_dist, n_references, norm, thresh, paired_end):
        self.clst_dist = clst_dist
        self.max_dist = max_dist
        self.local_chrom = -1
        self.norm = norm
        self.thresh = thresh
        self.paired_end = paired_end
        cdef cpp_map[int, LocalVal] scope
        for n in range(n_references + 1):  # Add one for special 'insertion chromosome'
            self.chrom_scope.push_back(scope)
    cdef void empty_scopes(self) nogil:
        for idx in range(self.chrom_scope.size()):
            if not self.chrom_scope[idx].empty():
                self.chrom_scope[idx].clear()
        self.loci.clear()

    cdef vector[int] find_other_nodes(self, int node_name, int current_chrom, int current_pos, int chrom2, int pos2,
                                      ReadEnum_t read_enum, int length_from_cigar, bint trust_ins_len) nogil:
        # todo make this code less nested and more readable
        cdef int idx, i, count_back, steps, node_name2
        cdef int sep = 0
        cdef int sep2 = 0
        cdef vector[int] found2
        cdef vector[int] found_exact
        cdef cpp_map[int, LocalVal]* forward_scope = &self.chrom_scope[chrom2]
        if chrom2 == 10000000:
            forward_scope = &self.chrom_scope.back()
        else:
            forward_scope = &self.chrom_scope[chrom2]
        cdef cpp_map[int, LocalVal].iterator local_it, local_it2
        cdef cpp_pair[int, LocalVal] vitem
        cdef float max_span, span_distance
        # Debug
        # echo("current_chrom ", current_chrom, "chrom2", chrom2, "current_pos", current_pos, "pos2", pos2)
        # echo("Forward scope len", forward_scope.size(), "Loci scope len", self.loci.size())
        # local_it = forward_scope.begin()
        # while local_it != forward_scope.end():
        #     vitem = dereference(local_it)
        #     echo(vitem.first, vitem.second)
        #     preincrement(local_it)

        # Re-initialize empty
        if current_chrom != self.local_chrom:
            self.local_chrom = current_chrom
            self.empty_scopes()
        if not self.loci.empty():
            # Erase items out of range in forward scope
            local_it = self.loci.lower_bound(current_pos - self.clst_dist)
            local_it2 = self.loci.begin()
            while local_it2 != local_it:
                vitem = dereference(local_it2)
                self.chrom_scope[vitem.second.chrom2].erase(vitem.second.pos2)
                preincrement(local_it2)
            if local_it != self.loci.begin():
                self.loci.erase(self.loci.begin(), local_it)

            local_it = forward_scope.lower_bound(pos2)
            steps = 0
            if local_it != forward_scope.end():
                while steps < 20: #6:
                    vitem = dereference(local_it)
                    if (read_enum == DELETION and vitem.second.read_enum == INSERTION) or (read_enum == INSERTION and vitem.second.read_enum == DELETION):
                        steps += 1
                        continue
                    node_name2 = vitem.second.node_name
                    if node_name2 != node_name:  # Can happen due to within-read events
                        # if node_name == 77:
                            # echo("-->", node_name, node_name2, current_chrom == chrom2, current_pos, pos2, pos2 - current_pos, vitem.second.pos2 - vitem.first,
                            #      is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2),
                            #      )
                            # echo(read_enum, vitem.second.read_enum, read_enum == INSERTION, vitem.second.read_enum == DELETION)
                        if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):
                            sep = c_abs(vitem.first - pos2)
                            sep2 = c_abs(vitem.second.pos2 - current_pos)
                            if sep < self.max_dist and sep2 < self.max_dist:
                                if sep < 35:
                                    if length_from_cigar > 0 and vitem.second.length_from_cigar > 0:
                                        max_span = max(length_from_cigar, vitem.second.length_from_cigar)
                                        span_distance = <float>c_abs(length_from_cigar - vitem.second.length_from_cigar) / max_span
                                        # echo("hi", node_name2, (length_from_cigar, vitem.second.length_from_cigar), span_distance, self.thresh)
                                        if span_distance < 0.8:
                                            found_exact.push_back(node_name2)
                                    else:
                                        found_exact.push_back(node_name2)
                                elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len):
                                    found2.push_back(node_name2)
                            elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len):
                                found2.push_back(node_name2)
                        elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len):
                            found2.push_back(node_name2)

                    if sep >= self.max_dist:
                        break
                    preincrement(local_it)
                    steps += 1
                    if local_it == forward_scope.end():
                        break

            if found_exact.empty():
                local_it = forward_scope.lower_bound(pos2)
                vitem = dereference(local_it)
                if local_it != forward_scope.begin():
                    predecrement(local_it)  # Move back one before staring search, otherwise same value is processed twice
                    steps = 0
                    while steps < 20: # 6:
                        vitem = dereference(local_it)
                        if (read_enum == DELETION and vitem.second.read_enum == INSERTION) or (read_enum == INSERTION and vitem.second.read_enum == DELETION):
                            steps += 1
                            continue
                        node_name2 = vitem.second.node_name
                        if node_name2 != node_name:
                            # if node_name == 1:
                            # echo("2", node_name, node_name2, read_enum) #, is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2), span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar))
                            # if node_name == 77:
                            #     echo('r', current_pos, pos2, vitem.first, vitem.second.pos2, span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len))
                            if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):
                                sep = c_abs(vitem.first - pos2)
                                sep2 = c_abs(vitem.second.pos2 - current_pos)
                                if sep < self.max_dist and vitem.second.chrom2 == chrom2 and \
                                        sep2 < self.max_dist:
                                    if sep < 35:
                                        if length_from_cigar > 0 and vitem.second.length_from_cigar > 0:
                                            max_span = max(length_from_cigar, vitem.second.length_from_cigar)
                                            span_distance = <float>c_abs(length_from_cigar - vitem.second.length_from_cigar) / max_span
                                            # echo("hi2", node_name2, (length_from_cigar, vitem.second.length_from_cigar), span_distance, self.thresh)
                                            if span_distance < 0.8:
                                                found_exact.push_back(node_name2)
                                        else:
                                            found_exact.push_back(node_name2)
                                    elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len):
                                        found2.push_back(node_name2)
                                elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len):
                                    found2.push_back(node_name2)
                            elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end, length_from_cigar, vitem.second.length_from_cigar, trust_ins_len):
                                found2.push_back(node_name2)
                        if local_it == forward_scope.begin() or sep >= self.max_dist:
                            break
                        predecrement(local_it)
                        steps += 1
        if not found_exact.empty():
            return found_exact
        else:
            return found2

    cdef void add_item(self, int node_name, int current_chrom, int current_pos, int chrom2, int pos2,
                       ReadEnum_t read_enum, int length_from_cigar): # nogil:

        # Add to scopes, if event is within read, add two references to forward scope. Otherwise when the
        # first break point drops from the local scope, the forward scope break will not connect with other
        # events that are added later. This is a fix for long read deletions mainly
        if chrom2 == -1:
            return  # No chrom2 was set, single-end?
        cdef cpp_map[int, LocalVal]* forward_scope
        if chrom2 == 10000000:
            forward_scope = &self.chrom_scope.back()
        else:
            forward_scope = &self.chrom_scope[chrom2]
        # Add to local scope
        cdef cpp_pair[int, LocalVal] pp
        pp.first = current_pos
        pp.second = make_local_val(chrom2, pos2, node_name, read_enum, length_from_cigar)
        self.loci.insert(pp)
        if read_enum == DELETION:
            forward_scope.insert(pp)
        # Add to forward scope
        pp.first = pos2
        pp.second = make_local_val(current_chrom, current_pos, node_name, read_enum, length_from_cigar)
        forward_scope.insert(pp)


cdef class TemplateEdges:
    cdef public unordered_map[string, vector[int]] templates_s
    def __init__(self):
        pass
    cdef void add(self, str template_name, int flag, int node, int query_start):
        cdef vector[int] val
        cdef bytes key = bytes(template_name, encoding="utf8")
        val.push_back(query_start)
        val.push_back(node)
        val.push_back(flag)
        self.templates_s[key].insert(self.templates_s[key].end(), val.begin(), val.end())


cdef void add_template_edges(Py_SimpleGraph G, TemplateEdges template_edges):
    # this function joins up template reads (read 1, read 2, plus any supplementary)
    cdef int ii, u_start, v_start, u, v, uflag, vflag
    cdef unordered_map[string, vector[int]].iterator it = template_edges.templates_s.begin()
    cdef vector[int] arr
    while it != template_edges.templates_s.end():
        arr = dereference(it).second  # Array values are query start, node-name, flag
        postincrement(it)
        read1_aligns = []
        read2_aligns = []
        for ii in range(0, arr.size(), 3):
            if arr[ii + 2] & 64:  # first in pair
                read1_aligns.append(arr[ii:ii + 3])
            else:
                read2_aligns.append(arr[ii:ii + 3])
        primary1 = None
        primary2 = None
        if len(read1_aligns) > 0:
            if len(read1_aligns) == 1:
                if not read1_aligns[0][2] & 2304:  # Is primary
                    primary1 = read1_aligns[0][1]
            else:
                if len(read1_aligns) > 2:
                    read1_aligns = sorted(read1_aligns)
                # Add edge between alignments that are neighbors on the query sequence, or between primary alignments
                for ii in range(len(read1_aligns) - 1):
                    u_start, u, uflag = read1_aligns[ii]
                    if not uflag & 2304:
                        primary1 = u
                    v_start, v, vflag = read1_aligns[ii + 1]
                    if not G.hasEdge(u, v):
                        G.addEdge(u, v, w=1)
                if primary1 is None:
                    if not read1_aligns[-1][2] & 2304:
                        primary1 = read1_aligns[-1][1]
        if len(read2_aligns) > 0:
            if len(read2_aligns) == 1:
                if not read2_aligns[0][2] & 2304:
                    primary2 = read2_aligns[0][1]
            else:
                if len(read2_aligns) > 2:
                    read2_aligns = sorted(read2_aligns)
                for ii in range(len(read2_aligns) - 1):
                    u_start, u, uflag = read2_aligns[ii]
                    if not uflag & 2304:
                        primary2 = u
                    v_start, v, vflag = read2_aligns[ii + 1]
                    if not G.hasEdge(u, v):
                        G.addEdge(u, v, w=1)
                if primary2 is None:
                    if not read2_aligns[-1][2] & 2304:
                        primary2 = read2_aligns[-1][1]
        if primary1 is not None and primary2 is not None:
            if not G.hasEdge(primary1, primary2):
                G.addEdge(primary1, primary2, w=1)


@cython.auto_pickle(True)
cdef class NodeName:
    cdef public uint64_t hash_name
    cdef public uint64_t tell
    cdef public uint32_t pos
    cdef public int32_t cigar_index
    cdef public uint32_t event_pos
    cdef public uint16_t flag
    cdef public uint16_t chrom
    def __init__(self, h, f, p, c, t, cigar_index, event_pos):
        self.hash_name = h
        self.flag = f
        self.pos = p
        self.chrom = c
        self.tell = t
        self.cigar_index = cigar_index
        self.event_pos = event_pos
    def as_tuple(self):
        return self.hash_name, self.flag, self.pos, self.chrom, self.tell, self.cigar_index, self.event_pos

cdef class NodeToName:
    # Index these vectors to get the unique 'template_name'
    # node names have the form (hash qname, flag, pos, chrom, tell, cigar index, event pos)
    cdef vector[uint64_t] h
    cdef vector[uint16_t] f
    cdef vector[uint32_t] p
    cdef vector[uint16_t] c
    cdef vector[uint64_t] t
    cdef vector[int32_t] cigar_index
    cdef vector[uint32_t] event_pos
    def __cinit__(self):
        pass
    cdef void append(self, long a, int b, int c, int d, long e, int f, int g) nogil:
        self.h.push_back(a)
        self.f.push_back(b)
        self.p.push_back(c)
        self.c.push_back(d)
        if e == -1:  # means to look in read buffer instead
            e = 0
        self.t.push_back(e)
        self.cigar_index.push_back(f)
        self.event_pos.push_back(g)
    def __getitem__(self, idx):
        return NodeName(self.h[idx], self.f[idx], self.p[idx], self.c[idx], self.t[idx], self.cigar_index[idx], self.event_pos[idx])
    cdef bint same_template(self, int query_node, int target_node):
        if self.h[query_node] == self.h[target_node]:
            return True

cdef get_query_pos_from_cigarstring(cigar, pos):
    # Infer the position on the query sequence of the alignment using cigar string
    cdef int end = 0
    cdef int start = 0
    cdef bint i = 0
    cdef int ref_end = pos
    cdef int slen
    for slen, opp in cigar:
        if not i and opp in "SH":
            start += slen
            end += slen
            i = 1
        elif opp == "D":
            ref_end += slen
        elif opp == "I":
            end += slen
        elif opp in "M=X":
            end += slen
            ref_end += slen
        i = 1
    return start, end, pos, ref_end


def get_query_pos_from_cigartuples(r):
    # Infer the position on the query sequence of the alignment using cigar string
    start = 0
    query_length = r.infer_read_length()  # Note, this also counts hard-clips
    end = query_length

    if r.cigartuples[0][0] == 4 or r.cigartuples[0][0] == 5:
        start += r.cigartuples[0][1]
    if r.cigartuples[-1][0] == 4 or r.cigartuples[-1][0] == 5:
        end -= r.cigartuples[-1][1]
    return start, end, query_length


AlnBlock = namedtuple("SA", ["query_start", "query_end", "ref_start", "ref_end", "chrom", "mq", "strand", "this"])
JoinEvent = namedtuple("JE", ["chrom", "event_pos", "chrom2", "pos2", "read_enum", "cigar_index"])


class AlignmentsSA:
    def __init__(self, r, gettid, paired_end=False):
        self.paired = paired_end
        self.index = None
        self.aln_strand = None
        self.query_aligns = []
        self.join_result = []
        self._alignments_from_sa(r, gettid)

    def connect_alignments(self, r, max_dist=1000, mq_thresh=0, read_enum=0):
        if len(self.query_aligns) > 1 and self.index is not None:
            if self.index > 0:
                self._connect_left(r, max_dist, mq_thresh, read_enum)
            if self.index < len(self.query_aligns) - 1:
                self._connect_right(r, max_dist, mq_thresh, read_enum)

    def _alignments_from_sa(self, r, gettid):
        qstart, qend, query_length = get_query_pos_from_cigartuples(r)
        this_aln = AlnBlock(query_start=qstart, query_end=qend,
                            ref_start=r.pos, ref_end=r.reference_end, chrom=r.rname, mq=r.mapq,
                            strand="-" if r.flag & 16 else "+", this=True)
        self.aln_strand = this_aln.strand
        query_aligns = [this_aln]
        for sa_block in r.get_tag("SA").split(";"):
            if sa_block == "":
                break
            sa = sa_block.split(",", 5)
            matches = [(int(slen), opp) for slen, opp in re.findall(r'(\d+)([A-Z]{1})', sa[3])]
            query_start, query_end, ref_start, ref_end = get_query_pos_from_cigarstring(matches, int(sa[1]))
            if this_aln.strand != sa[2]:
                start_temp = query_length - query_end
                query_end = start_temp + query_end - query_start
                query_start = start_temp
            query_aligns.append(AlnBlock(query_start, query_end, ref_start, ref_end, gettid(sa[0]), int(sa[4]), sa[2], False))
        query_aligns = sorted(query_aligns)
        cdef int index = 0
        for i, item in enumerate(query_aligns):
            if item.this:
                self.index = i
                break
        self.query_aligns = query_aligns

    def _connect_left(self, r, max_dist, mq_thresh, read_enum): #, max_gap_size=100):
        a = self.query_aligns[self.index]
        b = self.query_aligns[self.index - 1]
        # gap = abs(a.query_start - b.query_end)
        # if gap > max_gap_size:
        #     return
        event_pos = a.ref_start
        chrom = a.chrom
        cigar_index = 0
        if b.strand == self.aln_strand:
            pos2 = b.ref_end
        else:
            pos2 = b.ref_start
        chrom2 = b.chrom
        if b.mq < mq_thresh:
            self.join_result.append(JoinEvent(chrom, event_pos, chrom, event_pos, ReadEnum_t.BREAKEND, cigar_index))
            return
        elif self.paired:
            if not r.flag & 2304 and r.flag & 2 and r.pnext <= event_pos and (chrom != chrom2 or abs(pos2 - event_pos) > max_dist):
                self.join_result.append(JoinEvent(chrom, event_pos, chrom, event_pos, ReadEnum_t.BREAKEND, cigar_index))
                return
        self.join_result.append(JoinEvent(chrom, event_pos, chrom2, pos2, read_enum, cigar_index))

    def _connect_right(self, r, max_dist, mq_thresh, read_enum): #, max_gap_size=100):
        a = self.query_aligns[self.index]
        b = self.query_aligns[self.index + 1]
        # gap = abs(b.query_start - a.query_end)
        # if gap > max_gap_size:
        #     return
        event_pos = a.ref_end
        chrom = a.chrom
        cigar_index = len(r.cigartuples) - 1
        if b.strand == self.aln_strand:
            pos2 = b.ref_start
        else:
            pos2 = b.ref_end
        chrom2 = b.chrom
        if b.mq < mq_thresh:
            self.join_result.append(JoinEvent(chrom, event_pos, chrom, event_pos, ReadEnum_t.BREAKEND, cigar_index))
            return
        elif self.paired:
            # If paired, and SA block fits between two normal primary alignments, and block is not local, then
            # ignore block and try and call insertion
            if not r.flag & 2304 and r.flag & 2 and r.pnext >= event_pos and (chrom != chrom2 or abs(pos2 - event_pos) > max_dist):  # is primary, is proper pair. Use primary mate info, not SA
                self.join_result.append(JoinEvent(chrom, event_pos, chrom, event_pos, ReadEnum_t.BREAKEND, cigar_index))
                return
        self.join_result.append(JoinEvent(chrom, event_pos, chrom2, pos2, read_enum, cigar_index))


cdef int cluster_clipped(Py_SimpleGraph G, r, ClipScoper_t clip_scope, chrom, pos, node_name):
    cdef int other_node
    cdef int count = 0
    cdef unordered_set[int] clustered_nodes
    clip_scope.update(r, node_name, chrom, pos, clustered_nodes)
    if not clustered_nodes.empty():
        for other_node in clustered_nodes:
            if not G.hasEdge(node_name, other_node):
                G.addEdge(node_name, other_node, 2)
                count += 1
    return count


cdef void add_to_graph(Py_SimpleGraph G, AlignedSegment r, PairedEndScoper_t pe_scope, TemplateEdges_t template_edges,
                       NodeToName node_to_name, genome_scanner,
                       int flag, int chrom, tell, int cigar_index, int event_pos,
                       int chrom2, int pos2, ClipScoper_t clip_scope, ReadEnum_t read_enum,
                       bint p1_overlaps, bint p2_overlaps, bint mm_only, int clip_l, site_adder,
                       int length_from_cigar, bint trust_ins_len, bint paired_end):
    # Adds relevant information to graph and other data structures for further processing
    cdef int other_node
    cdef vector[int] other_nodes  # Other alignments to add edges between
    cdef int node_name = G.addNode()
    cdef uint64_t v = xxhasher(bam_get_qname(r._delegate), len(r.qname), 42)  # Hash qname to save mem
    cdef int bnd_site_node, bnd_site_node2
    cdef int q_start
    # echo("\nADDING", node_name, r.qname, (chrom, event_pos), (chrom2, pos2), flag, read_enum)
    node_to_name.append(v, flag, r.pos, chrom, tell, cigar_index, event_pos)  # Index this list to get the template_name
    genome_scanner.add_to_buffer(r, node_name, tell)  # Add read to buffer
    if read_enum < 2:  # Prevents joining up within-read svs with between-read svs
        q_start = r.query_alignment_start if not r.flag & 16 else r.infer_query_length() - r.query_alignment_end
        template_edges.add(r.qname, flag, node_name, q_start)

    both_overlap = p1_overlaps and p2_overlaps
    if not paired_end or (paired_end and read_enum != BREAKEND and not mm_only and not chrom2 == -1):
        other_nodes = pe_scope.find_other_nodes(node_name, chrom, event_pos, chrom2, pos2, read_enum, length_from_cigar, trust_ins_len)
        for other_node in other_nodes:
            if node_to_name.same_template(node_name, other_node):
                continue
            if not G.hasEdge(node_name, other_node):
                G.addEdge(node_name, other_node, 2)  # 'black' edge

    elif chrom != chrom2 and clip_l != -1:  # Note all paired-end reads have BREAKENDS where chrom != chrom2, but also includes translocations
        cluster_clipped(G, r, clip_scope, chrom, event_pos, node_name)

    if not paired_end or (paired_end and not mm_only and read_enum != BREAKEND):
        pe_scope.add_item(node_name, chrom, event_pos, chrom2, pos2, read_enum, length_from_cigar)
    if site_adder and (read_enum == BREAKEND or read_enum == SPLIT):
        bnd_site_node = site_adder.find_nearest_site(chrom, event_pos)
        if bnd_site_node >= 0 and not G.hasEdge(node_name, bnd_site_node):
            G.addEdge(node_name, bnd_site_node, 0)
        bnd_site_node2 = site_adder.find_nearest_site(chrom2, pos2)
        if bnd_site_node != bnd_site_node2 and bnd_site_node2 >= 0 and not G.hasEdge(node_name, bnd_site_node2):
            G.addEdge(node_name, bnd_site_node2, 0)
    # Debug:
    # if r.qname == "M03762:232:000000000-L65J4:1:2111:17729:15161":
    #     echo("---", r.qname, read_enum, node_name, (event_pos, pos2), length_from_cigar, list(other_nodes))
    # look = {'a4d38568-fd80-4785-8fa5-84ed132b445c', '2313a985-385c-4c84-b02c-dddfc627940b', '0031840a-bd2d-475d-9a04-528f71c7b512'}
    # if r.qname in look:
    # # if r.qname == "D00360:18:H8VC6ADXX:1:1210:7039:44052":
    #     echo(r.qname, r.flag, node_name, (chrom, event_pos), (chrom2, pos2), list(other_nodes),
    #          cigar_index, length_from_cigar, "enum=", read_enum)
    # echo(list(other_nodes))


cdef int good_quality_clip(AlignedSegment r, int clip_length):
    # Use sliding window to check that soft-clip has at least n good bases over qual 20
    cdef const unsigned char[:] quals = r.query_qualities
    if len(quals) == 0:
        return 1
    cdef char* char_ptr_rseq = <char*>bam_get_seq(r._delegate)
    cdef int i, length, w_sum
    cdef int window_length = 10
    if clip_length < window_length:
        window_length = clip_length
    cdef float avg_qual
    cdef int poly = window_length - 1
    cdef total_good = window_length - 1
    ct = r.cigartuples
    if r.flag & 2304:  # supplementary, usually hard-clipped, assume good clipping
        if ct[0][0] == 5 or ct[-1][0] == 5:
            return 1
    cdef int[17] basecounts  # char value is index into array
    cdef bint homopolymer
    cdef int base_index
    if ct[0][0] == 4:
        length = r.cigartuples[0][1]
        if length >= window_length and length >= clip_length:
            for i in range(length - window_length, -1, -1):
                # average of current window by counting leftwards
                for j in range(0, 16):
                    basecounts[j] = 0
                homopolymer = False
                w_sum = 0
                for j in range(i, i + window_length):
                    w_sum += quals[j]
                    base_index = <int>bam_seqi(char_ptr_rseq, j)
                    basecounts[base_index] += 1
                    if basecounts[base_index] == poly:
                        homopolymer = True
                        break
                avg_qual = w_sum / float(window_length)
                if avg_qual > 10:
                    if not homopolymer:
                        total_good += 1
                else:
                    break
            if total_good >= clip_length:
                return 1

    total_good = window_length - 1
    if ct[-1][0] == 4:
        length = r.cigartuples[-1][1]
        if length >= window_length and length >= clip_length:
            for i in range(len(r.query_qualities) - length, len(r.query_qualities) - window_length):
                # average of current window by counting rightwards
                for j in range(0, 16):
                    basecounts[j] = 0
                homopolymer = False
                w_sum = 0
                for j in range(i, i + window_length):
                    w_sum += quals[j]
                    base_index = <int>bam_seqi(char_ptr_rseq, j)
                    basecounts[base_index] += 1
                    if basecounts[base_index] == poly:
                        homopolymer = True
                        break
                avg_qual = w_sum / float(window_length)
                if avg_qual > 10:
                    if not homopolymer:
                        total_good += 1
                else:
                    break
            if total_good >= clip_length:
                return 1
    return 0


cdef void process_alignment(Py_SimpleGraph G, AlignedSegment r, int clip_l, int loci_dist, gettid,
                            overlap_regions, int clustering_dist, PairedEndScoper_t pe_scope,
                            int cigar_index, int event_pos, int paired_end, long tell, genome_scanner,
                            TemplateEdges_t template_edges, NodeToName node_to_name,
                            int cigar_pos2, int mapq_thresh, ClipScoper_t clip_scope,
                            ReadEnum_t read_enum, bad_clip_counter, bint mm_only, site_adder,
                            int length_from_cigar, bint trust_ins_len):
    cdef int other_node, clip_left, clip_right
    cdef bint current_overlaps_roi, next_overlaps_roi
    cdef bint add_primark_link, add_insertion_link
    cdef int chrom = r.rname
    cdef int chrom2 = r.rnext
    cdef int pos2 = r.pnext
    cdef int flag = r.flag
    cdef int left_clip, right_clip
    cdef str qname = r.qname
    cdef ReadEnum_t inferred_clip_type
    cdef uint64_t v
    cdef bint success
    cdef bint good_clip
    if paired_end and read_enum == SPLIT and flag & 8:  # clip event, or whole read, but mate is unmapped
        return
    if site_adder:
        if site_adder:
            site_adder.add_any_sites(r.rname, event_pos, G, pe_scope, node_to_name, clustering_dist)
    if paired_end or flag & 1:
        pnext = r.pnext
        rnext = r.rnext
        chrom2 = rnext
        pos2 = pnext
        current_overlaps_roi = intersecter(overlap_regions, r.rname, event_pos, event_pos+1)
        next_overlaps_roi = False
        if current_overlaps_roi:  # check if both sides of SV are in overlaps regions
            next_overlaps_roi = intersecter(overlap_regions, chrom2, pos2, pnext+1)
        # Special treatment of supplementary and local reads; need to decide where the partner is
        # Either use the rnext:pnext or a site(s) listed in the SA tag: The rnext:pext can often be at the
        # same loci as the read which leads to problems when linking w=2 edges
        add_primary_link = 1
        add_insertion_link = 0
        good_clip = good_quality_clip(r, 15)
        if not good_clip:
            bad_clip_counter.add(chrom, r.pos)
        if read_enum == SPLIT:
            if r.has_tag("SA") and good_clip:
                all_aligns = AlignmentsSA(r, gettid, True)
                all_aligns.connect_alignments(r, loci_dist, mapq_thresh, read_enum)
                if not all_aligns.join_result:
                    return
                for j in all_aligns.join_result:
                    chrom = j.chrom
                    event_pos = j.event_pos
                    chrom2 = j.chrom2
                    pos2 = j.pos2
                    read_enum = j.read_enum
                    cigar_index = j.cigar_index
                    if read_enum == BREAKEND:
                        if mm_only and current_overlaps_roi and next_overlaps_roi:  # Probably too many reads in ROI to reliably separate out break end reads
                            return
                        if good_clip and cigar_clip(r, clip_l):
                            chrom2 = 10000000
                            pos2 = event_pos
                        else:
                            return
                    if read_enum == DELETION or read_enum == INSERTION:
                        chrom2 = chrom
                        if r.cigartuples[cigar_index][0] != 1:  # not insertion, use length of cigar event
                            pos2 = cigar_pos2
                        else:
                            pos2 = event_pos

                    add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                 tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum,
                                 current_overlaps_roi, next_overlaps_roi,
                                 mm_only, clip_l, site_adder, 0, trust_ins_len, paired_end)
                return

        if read_enum == BREAKEND:
            if mm_only and current_overlaps_roi and next_overlaps_roi:  # Probably too many reads in ROI to reliably separate out break end reads
                return
            if good_clip and cigar_clip(r, clip_l):
                chrom2 = 10000000
                pos2 = event_pos
            else:
                return
        if read_enum == DELETION or read_enum == INSERTION:
            chrom2 = chrom
            if r.cigartuples[cigar_index][0] != 1:  # not insertion, use length of cigar event
                pos2 = cigar_pos2
            else:
                pos2 = event_pos
        add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                               tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum, current_overlaps_roi, next_overlaps_roi,
                               mm_only, clip_l, site_adder, length_from_cigar, trust_ins_len, paired_end)
    ###
    else:  # Single end
        current_overlaps_roi, next_overlaps_roi = False, False  # not supported
        if read_enum == SPLIT:
            if r.has_tag("SA"):
                all_aligns = AlignmentsSA(r, gettid, False)
                all_aligns.connect_alignments(r, loci_dist, mapq_thresh, read_enum)
                if not all_aligns.join_result:
                    return
                # for chrom, event_pos, chrom2, pos2, read_e, cigar_index in all_aligns.join_result:
                for j in all_aligns.join_result:
                    add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, j.chrom,
                                 tell, j.cigar_index, j.event_pos, j.chrom2, j.pos2, clip_scope, j.read_enum,
                                 current_overlaps_roi, next_overlaps_roi,
                                 mm_only, clip_l, site_adder, 0, trust_ins_len, paired_end)

        elif read_enum >= 2:  # Sv within read or breakend
            chrom2 = r.rname
            if cigar_index != -1 and r.cigartuples[cigar_index][0] != 1:  # If not insertion
                pos2 = cigar_pos2
            else:
                pos2 = event_pos
            add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                         tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum,
                         current_overlaps_roi, next_overlaps_roi, mm_only, clip_l, site_adder, length_from_cigar, trust_ins_len, paired_end)


cdef struct CigarEvent:
    int opp
    int cigar_index
    int event_pos
    int pos2
    int length
    bint cigar_skip
    ReadEnum_t read_enum


cdef CigarEvent make_cigar_event(int opp, int cigar_index, int event_pos, int pos2, int length,
                                 ReadEnum_t read_enum) nogil:
    cdef CigarEvent item
    item.opp = opp
    item.cigar_index = cigar_index
    item.event_pos = event_pos
    item.read_enum = read_enum
    item.pos2 = pos2
    item.length = length
    item.cigar_skip = 0
    return item


class SiteAdder:
    def __init__(self, sites):
        self.sites = {k: list(v) for k, v in sites.items()}
        self.sites_queue = sites  # dict of deque's
        self.sites_scope = sites
        self.sites_index = {}
        self.sites_index_r = {}
        self.current_index = 0
        self.count = 0
        self.scope = sortedcontainers.sortedlist.SortedList(key=lambda x: x.start)
        self.current_chrom = -1
        self.scope_front = 0
        self.scope_back = 0
        self.lookup = namedtuple("l", ["start"])
    def __contains__(self, item):
        if item in self.sites_index:
            return True
    def __getitem__(self, item):
        return self.sites_index[item]
    def find_nearest_site(self, chrom, pos):
        # search backwards from current_index for closest --site to pos
        cluster_dist = 50
        scope_dist = 500
        if len(self.scope) == 0:
            return -1
        left_index = self.scope.bisect_left(self.lookup(pos - scope_dist)) - 1
        if left_index > 0:
            for i in range(left_index):
                self.scope.pop(0)
        left_index = self.scope.bisect_left(self.lookup(pos))
        if left_index == len(self.scope):
            p = self.scope[left_index - 1]
        elif left_index != len(self.scope) - 1:
            left_p = self.scope[left_index]
            right_p = self.scope[left_index + 1]
            left_dis = abs(left_p.start - pos)
            right_dis = abs(right_p.start - pos)
            if left_dis <= right_dis:
                p = left_p
            else:
                p = right_p
        else:
            p = self.scope[left_index]
        if abs(p.start - pos) <= cluster_dist:
            return self.sites_index_r[p]
        else:
            return -1
    def add_to_scope(self, site):
        if site.chrom == self.current_chrom:
            self.scope.add(site)
        else:
            self.scope = sortedcontainers.sortedlist.SortedList([site], key=lambda x: x.start)
            self.current_chrom = site.chrom
    def add_any_sites(self, int chrom, int pos, Py_SimpleGraph G, PairedEndScoper_t pe_scope, NodeToName node_to_name, cluster_dist):
        cdef int node_name, start, stop, file_index
        cdef ReadEnum_t read_enum
        if chrom not in self.sites_queue:
            return
        s = self.sites_queue[chrom]
        if len(s) == 0:
            return
        p = s[0]
        if pos + cluster_dist < p.start:
            return
        while p.start < pos - cluster_dist and len(s) > 0:
            s.popleft()
            if len(s) == 0:
                return
            p = s[0]
        while abs(p.start - pos) < cluster_dist and len(s) > 0:
            site = p
            if site.svtype == "DEL":
                read_enum = DELETION
                length = site.svlen
            elif site.svtype == "INS":
                read_enum = INSERTION
                length = site.svlen
            else:
                read_enum = BREAKEND
                length = 0  # cluster by distance only
            node_name = G.addNode()
            self.sites_index[node_name] = site
            self.sites_index_r[site] = node_name
            pe_scope.add_item(node_name, site.chrom, site.start, site.chrom2, site.end, read_enum, length)
            pe_scope.local_chrom = chrom  # needs setting otherwise scope can be cleared when new read is added
            node_to_name.append(0, 0, 0, 0, 0, 0, 0)
            self.count += 1
            self.sites_queue[chrom].popleft()
            self.add_to_scope(site)
            if len(s) > 0:
                p = s[0]
            else:
                break


cpdef tuple construct_graph(genome_scanner, infile, int max_dist, int clustering_dist, int k=16, int m=7, int clip_l=21,
                            int min_sv_size=30,
                            int minimizer_support_thresh=2, int minimizer_breadth=3,
                            int minimizer_dist=10, int mapq_thresh=1, debug=None, procs=1,
                            int paired_end=1, int read_length=150, bint contigs=True,
                            float norm_thresh=100, float spd_thresh=0.3, bint mm_only=False,
                            sites=None, bint trust_ins_len=True, low_mem=False, temp_dir=".",
                            find_n_aligned_bases=True):
    logging.info("Building graph with clustering {} bp".format(clustering_dist))
    cdef TemplateEdges_t template_edges = TemplateEdges()  # Edges are added between alignments from same template, after building main graph
    cdef int event_pos, cigar_index, opp, length
    node_to_name = NodeToName()  # Map of nodes -> read ids
    cdef ClipScoper_t clip_scope = ClipScoper(minimizer_dist, k=k, m=m, clip_length=clip_l,  # Keeps track of local reads
                       minimizer_support_thresh=minimizer_support_thresh,
                       minimizer_breadth=minimizer_breadth, read_length=read_length)
    # Infers long-range connections, outside local scope using pe information
    cdef PairedEndScoper_t pe_scope = PairedEndScoper(max_dist, clustering_dist, infile.header.nreferences, norm_thresh, spd_thresh, paired_end)
    bad_clip_counter = BadClipCounter(infile.header.nreferences, low_mem, temp_dir)
    cdef Py_SimpleGraph G = map_set_utils.Py_SimpleGraph()
    site_adder = None
    if sites:
        site_adder = SiteAdder(sites)
    overlap_regions = genome_scanner.overlap_regions  # Get overlapper, intersect reads with intervals
    gettid = infile.gettid
    cdef long tell
    cdef int pos2
    cdef bint added
    cdef ReadEnum_t read_enum
    cdef bint clipped
    cdef vector[CigarEvent] events_to_add
    cdef vector[CigarEvent].iterator itr_events
    cdef CigarEvent v
    cdef AlignedSegment r
    cdef uint32_t cigar_value
    cdef uint32_t cigar_l
    cdef uint32_t *cigar_p
    cdef long n_aligned_bases = 0

    for chunk in genome_scanner.iter_genome():
        for r, tell in chunk:
            if r.mapq < mapq_thresh:
                continue
            pos2 = -1
            event_pos = r.pos
            added = False
            clipped = 0
            events_to_add.clear()
            cigar_l = r._delegate.core.n_cigar
            cigar_p = bam_get_cigar(r._delegate)
            if cigar_l > 1:
                if r.has_tag("SA"):
                    # Set cigar-index to -1 means it is unset, will be determined during SA parsing
                    cigar_index = -1
                    process_alignment(G, r, clip_l, max_dist, gettid,
                                      overlap_regions, clustering_dist, pe_scope,
                                      cigar_index, event_pos, paired_end, tell, genome_scanner,
                                      template_edges, node_to_name, pos2, mapq_thresh, clip_scope, ReadEnum_t.SPLIT,
                                      bad_clip_counter, mm_only, site_adder, 0, trust_ins_len)
                    added = True
                for cigar_index in range(cigar_l):
                    cigar_value = cigar_p[cigar_index]
                    opp = <int> cigar_value & 15
                    length = <int> cigar_value >> 4

                    if find_n_aligned_bases and (opp == 0 or opp == 7 or opp == 8):
                        n_aligned_bases += length

                    if opp == 1:
                        if length >= min_sv_size:
                            # Insertion type
                            pos2 = event_pos + length
                            events_to_add.push_back(make_cigar_event(opp, cigar_index, event_pos, pos2, length, ReadEnum_t.INSERTION))
                            added = True
                    elif opp == 2:
                        if length >= min_sv_size:  # elif!
                            pos2 = event_pos + length
                            events_to_add.push_back(make_cigar_event(opp, cigar_index, event_pos, pos2, length, ReadEnum_t.DELETION))
                            added = True
                        event_pos += length
                    else:
                        if opp != 4 and opp != 5:
                            event_pos += length
                        elif opp == 4 and length >= clip_l:
                            clipped = 1

            if (paired_end and not added and contigs) or not paired_end:
                # Whole alignment will be used, try infer position from soft-clip
                cigar_index = -1
                pos2 = -1
                left_clip_size, right_clip_size = clip_sizes_hard(r)  # soft and hard-clips
                if r.flag & 8 and clipped:  # paired by inference
                    # skip if both ends are clipped, usually means its a chunk of badly mapped sequence
                    # if not (left_clip_size and right_clip_size) and ((paired_end and good_quality_clip(r, 20)) or (not paired_end and) ):
                    if not (left_clip_size and right_clip_size) and good_quality_clip(r, 20):
                        # Mate is unmapped, insertion type. Only add if soft-clip is available
                        process_alignment(G, r, clip_l, max_dist, gettid,
                                          overlap_regions, clustering_dist, pe_scope,
                                          cigar_index, event_pos, paired_end, tell, genome_scanner,
                                          template_edges, node_to_name,
                                          pos2, mapq_thresh, clip_scope, ReadEnum_t.BREAKEND, bad_clip_counter,
                                          mm_only, site_adder, 0, trust_ins_len)
                else:
                    # Use whole read, could be normal or discordant
                    if not paired_end:
                        if max(left_clip_size, right_clip_size) > 250:
                            read_enum = ReadEnum_t.BREAKEND
                            if left_clip_size > right_clip_size:
                                event_pos = r.pos  # else reference_end is used
                            process_alignment(G, r, clip_l, max_dist, gettid,
                                              overlap_regions, clustering_dist, pe_scope,
                                              cigar_index, event_pos, paired_end, tell, genome_scanner,
                                              template_edges, node_to_name,
                                              pos2, mapq_thresh, clip_scope, read_enum, bad_clip_counter,
                                              mm_only, site_adder, 0, trust_ins_len)

                    else:
                        if r.flag & 2 and abs(r.tlen) < max_dist and r.rname == r.rnext:
                            if not clipped:
                                continue
                            read_enum = ReadEnum_t.BREAKEND
                        else:
                            read_enum = ReadEnum_t.DISCORDANT

                        if left_clip_size or right_clip_size:
                            if left_clip_size > right_clip_size:
                                event_pos = r.pos  # else reference_end is used
                        process_alignment(G, r, clip_l, max_dist, gettid,
                                          overlap_regions, clustering_dist, pe_scope,
                                          cigar_index, event_pos, paired_end, tell, genome_scanner,
                                          template_edges, node_to_name,
                                          pos2, mapq_thresh, clip_scope, read_enum, bad_clip_counter,
                                          mm_only, site_adder, 0, trust_ins_len)
            # process within-read events
            if not events_to_add.empty():
                itr_events = events_to_add.begin()
                while itr_events != events_to_add.end():
                    v = dereference(itr_events)
                    if v.cigar_skip:
                        pos2 = v.pos2
                    else:
                        pos2 = v.event_pos + v.length  # fall back on original cigar event length
                    process_alignment(G, r, clip_l, max_dist, gettid,
                                      overlap_regions, clustering_dist, pe_scope,
                                      v.cigar_index, v.event_pos, paired_end, tell, genome_scanner,
                                      template_edges, node_to_name,
                                      v.pos2, mapq_thresh, clip_scope, v.read_enum, bad_clip_counter,
                                      mm_only, site_adder, v.length, trust_ins_len)
                    preincrement(itr_events)

    add_template_edges(G, template_edges)
    if site_adder:
        logging.info(f"Added {site_adder.count} variants from input sites")
    return G, node_to_name, bad_clip_counter, site_adder, n_aligned_bases


cdef BFS_local(Py_SimpleGraph G, int source, unordered_set[int]& visited ):
    cdef array.array queue = array.array("L", [source])
    nodes_found = set([])
    cdef int u, v
    cdef vector[int] neighbors
    while queue:
        u = queue.pop(0)
        G.neighbors(u, neighbors)
        for v in neighbors:
            if visited.find(v) == visited.end():
                if G.weight(u, v) > 1:
                    if u not in nodes_found:
                        nodes_found.add(u)
                    if v not in nodes_found:
                        nodes_found.add(v)
                        queue.append(v)
        visited.insert(u)
    return array.array("L", nodes_found)


cdef get_partitions(Py_SimpleGraph G, nodes):
    cdef unordered_set[int] seen
    cdef int u, v, i
    cdef vector[int] neighbors
    parts = []
    for u in nodes:
        if seen.find(u) != seen.end():
            continue
        G.neighbors(u, neighbors)
        for v in neighbors:
            if seen.find(v) != seen.end():
                continue
            if G.weight(u, v) > 1:  # weight 2 or 3 for normal or black edges
                found = BFS_local(G, u, seen)
                if len(found):
                    parts.append(found)
        seen.insert(u)
    return parts


cdef tuple count_support_between(Py_SimpleGraph G, parts):
    cdef int i, j, node, child, any_out_edges
    cdef tuple t
    cdef unsigned long[:] p
    if len(parts) == 0:
        return {}, {}
    elif len(parts) == 1:
        return {}, {0: parts[0]}
    cdef Py_Int2IntMap p2i = map_set_utils.Py_Int2IntMap()
    for i, p in enumerate(parts):
        for node in p:
            p2i.insert(node, i)
    # Count the links between partitions. Split reads into sets ready for calling
    # No counting of read-pairs templates or 'support', just a count of linking alignments
    # counts (part_a, part_b): {part_a: {node 1, node 2 ..}, part_b: {node4, ..} }
    counts = {}
    self_counts = {}
    seen_t = set([])
    cdef vector[int] neighbors
    for i, p in enumerate(parts):
        current_t = set([])
        for node in p:
            any_out_edges = 0
            G.neighbors(node, neighbors)
            for child in neighbors:
                if not p2i.has_key(child):
                    continue  # Exterior child, not in any partition
                j = p2i.get(child)  # Partition of neighbor node
                if j != i:
                    any_out_edges = 1
                    if j < i:
                        t = (j, i)
                    else:
                        t = (i, j)
                    if t in seen_t:
                        continue
                    if t not in counts:
                        counts[t] = [set([]), set([])]
                    if j < i:
                        counts[t][0].add(child)
                        counts[t][1].add(node)
                    else:
                        counts[t][1].add(child)
                        counts[t][0].add(node)
                    current_t.add(t)
            # Count self links, important for resolving small SVs
            if any_out_edges == 0:
                if i not in self_counts:
                    self_counts[i] = array.array("L", [node])
                else:
                    self_counts[i].append(node)
        seen_t.update(current_t)  # Only count edge once
        for t in current_t:  # save memory by converting support_between to array
            counts[t] = [np.fromiter(m, dtype="uint32", count=len(m)) for m in counts[t]]
    return counts, self_counts


cpdef break_large_component(Py_SimpleGraph G, component, int min_support):
    # similar to count_support_between except only counts are returned without the node labels partitioned
    parts = get_partitions(G, component)
    cdef int i, j, node, child
    cdef tuple t
    if len(parts) <= 1:
        return parts
    # Make a table to count from, int-int
    cdef Py_Int2IntMap p2i = map_set_utils.Py_Int2IntMap()
    for i, p in enumerate(parts):
        for node in p:
            p2i.insert(node, i)
    # Count the links between partitions. Split reads into sets ready for calling
    # No counting of read-pairs templates or 'support', just a count of linking alignments
    # counts (part_a, part_b): number-of-links
    counts = defaultdict(int)
    self_counts = defaultdict(int)
    seen_t = set([])
    cdef vector[int] neighbors
    for i, p in enumerate(parts):
        current_t = set([])
        for node in p:
            any_out_edges = 0  # Keeps track of number of outgoing pairs, or self edges
            G.neighbors(node, neighbors)
            for child in neighbors:
                if not p2i.has_key(child):
                    continue  # Exterior child, not in any partition
                # Partition of neighbor node
                j = p2i.get(child)
                if j != i:
                    if j < i:
                        t = (j, i)
                    else:
                        t = (i, j)
                    if t in seen_t:
                        continue
                    counts[t] += 1
                    current_t.add(t)
                else:
                    self_counts[i] += 1
        seen_t.update(current_t)  # Only count edge once
    f = set([])
    jobs = []
    for (u, v), sup in counts.items():
        if sup >= min_support:
            f.add(u)
            f.add(v)
            if len(parts[u]) > 500 or len(parts[v]) > 500:
                continue
            b = list(parts[u])
            b.extend(parts[v])
            jobs.append(b)
    for k, v in self_counts.items():
        if v >= min_support and k not in f:
            jobs.append(parts[k])
    return jobs


cpdef proc_component(node_to_name, component, read_buffer, infile, Py_SimpleGraph G, int min_support, int procs, int paired_end,
                     sites_index):
    n2n = {}
    reads = {}
    cdef int support_estimate = 0
    cdef int v
    info = None
    if min_support >= 3 and len(component) == 1 and not sites_index:
        return
    for v in component:
        # Add information from --sites, keep any read found that link to --sites variants
        if sites_index and v in sites_index:
            if info is None:
                info = {v: sites_index[v]}
            else:
                info[v] = sites_index[v]
            min_support = len(info) + 1
            continue
        if procs == 1 and v in read_buffer:
            reads[v] = read_buffer[v]
        key = node_to_name[v]
        if key.cigar_index != -1:
            support_estimate += 2
        else:
            support_estimate += 1
        n2n[v] = key
    if support_estimate < min_support:
        return

    # Explore component for locally interacting nodes; create partitions using these
    partitions = get_partitions(G, component)
    support_between, support_within = count_support_between(G, partitions)
    if len(support_between) == 0 and len(support_within) == 0:
        if not paired_end:
            if len(n2n) >= min_support or len(reads) >= min_support or info:
                d = {"parts": [], "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
                if info:
                    d["info"] = info
                return d
            else:
                return
        else:
            # single paired end template can have 3 nodes e.g. two reads plus supplementary
            if min_support == 1 and (len(n2n) >= min_support or len(reads) >= min_support):
                d = {"parts": [], "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
                if info:
                    d["info"] = info
                return d
            elif len(reads) >= min_support or info:
                d = {"parts": [], "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
                if info:
                    d["info"] = info
                return d
            else:
                return
    # Debug:
    # if 157835 in n2n:
    #     echo("parts", partitions)
    #     echo("s_between", support_between)
    #     echo("s_within", support_within)

    # echo("parts", partitions)
    # echo("s_between", support_between)
    # echo("s_within", support_within)

    # echo("n2n", n2n.keys())
    # node_look = set(range(653526, 653532))
    # node_look = set(range(8))
    # echo(node_look.intersection(set(n2n.keys())))
    # if node_look.intersection(set(n2n.keys())):
    #     echo("parts", partitions)
    #     echo("s_between", sb)
    #     echo("s_within", support_within)
    d = {"parts": partitions, "s_between": support_between, "reads": reads, "s_within": support_within, "n2n": n2n}
    if info:
        d["info"] = info
    return d
