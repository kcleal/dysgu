#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

from __future__ import absolute_import
import time
import click
import ncls
from collections import defaultdict, deque
import numpy as np
cimport numpy as np
from cpython cimport array
import array
import re
import logging
import dataclasses

from libcpp.vector cimport vector
from libcpp.pair cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.unordered_map cimport unordered_map
from dysgu.map_set_utils cimport unordered_map as robin_map

from libcpp.string cimport string
from libcpp.deque cimport deque as cpp_deque

from libc.stdint cimport uint16_t, uint32_t, int32_t, uint64_t
from libc.stdlib cimport abs as c_abs

from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement

from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_qname, bam_seqi, bam_get_seq

from dysgu cimport map_set_utils
from dysgu import io_funcs
from dysgu.map_set_utils cimport unordered_set, cigar_clip, clip_sizes_hard, is_reciprocal_overlapping, span_position_distance, position_distance
from dysgu.map_set_utils cimport hash as xxhasher
from dysgu.map_set_utils cimport MinimizerTable
from dysgu.post_call_metrics import BadClipCounter

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


def echo(*args):
    click.echo(args, err=True)


cdef class Table:
    # overlap table for ncls
    cdef vector[np.int64_t] starts
    cdef vector[np.int64_t] ends
    cdef vector[np.int64_t] values

    cpdef void add(self, int s, int e, int v):
        self.starts.push_back(s)
        self.ends.push_back(e)
        self.values.push_back(v)

    def get_val(self, v):
        cdef vector[np.int64_t] values = v
        cdef np.ndarray[np.int64_t] a = np.empty(values.size(), dtype=np.int)
        cdef int len_a = len(a)
        cdef int i
        with nogil:
            for i in range(len_a):
                a[i] = values[i]
        return a

    def containment_list(self):
        return ncls.NCLS(self.get_val(self.starts), self.get_val(self.ends), self.get_val(self.values))


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
    cdef int n_local_minimizers
    cdef object rds_clip
    cdef bint minimizer_mode

    """Keeps track of which reads are in scope"""
    def __cinit__(self, int max_dist, int k, int m, int clip_length, int minimizer_support_thresh,
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

        self.n_local_minimizers = 0
        self.target_density = 2. / (m + 1)
        self.upper_bound_n_minimizers = read_length * self.target_density

        self.rds_clip = [deque([]), deque([])]

        self.minimizer_mode = False

    cdef void _add_m_find_candidates(self, clip_seq, int name, int idx, int position,
                                     robin_map[uint64_t, unordered_set[uint64_t]]& clip_table,
                                     # clip_table,
                                     unordered_set[int]& clustered_nodes):

        cdef unordered_set[long] clip_minimizers
        cdef unordered_set[uint64_t].iterator targets_iter, targets_end

        sliding_window_minimum(self.k, self.w, clip_seq, clip_minimizers) #, self.homopolymer_table)  # get minimizers of sequence
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

        n_local_minimizers = self.n_local_minimizers
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
                self.n_local_minimizers += 1
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
                             robin_map[uint64_t, unordered_set[uint64_t]]& clip_table,
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
                                self.n_local_minimizers -= 1

                        preincrement(set_iter)
                    mm_table.erase(name)
            else:
                break

    cdef void _insert(self, seq, int cigar_start, int cigar_end, int input_read, int position,
                      unordered_set[int]& clustered_nodes):
        # Find soft-clips of interest
        if cigar_start >= self.clip_length:
            clip_seq = left_soft_clips(seq, cigar_start)
            self._refresh_scope(self.scope_left, position, self.read_minimizers_left, self.clip_table_left)
            self._add_m_find_candidates(clip_seq, input_read, 0, position, self.clip_table_left, clustered_nodes)
            self.scope_left.push_back(cpp_item(position, input_read))

        if cigar_end >= self.clip_length:
            clip_seq = right_soft_clips(seq, cigar_end)
            self._refresh_scope(self.scope_right, position, self.read_minimizers_right, self.clip_table_right)
            self._add_m_find_candidates(clip_seq, input_read, 1, position, self.clip_table_right, clustered_nodes)
            self.scope_right.push_back(cpp_item(position, input_read))

    cdef void update(self, r, int input_read, int chrom, int position,
               unordered_set[int]& clustered_nodes):

        clip_left, clip_right = map_set_utils.clip_sizes(r)
        # Find candidate nodes with matching minimizers
        if chrom != self.current_chrom:
            self.scope_left.clear()
            self.scope_right.clear()
            self.clip_table_left.clear()
            self.clip_table_right.clear()
            self.current_chrom = chrom
        self._insert(r.seq, clip_left, clip_right, input_read, position, clustered_nodes)


cdef struct LocalVal:
    int chrom2
    int pos2
    int node_name
    ReadEnum_t read_enum


cdef LocalVal make_local_val(int chrom2, int pos2, int node_name, ReadEnum_t read_enum) nogil:
    cdef LocalVal item
    item.chrom2 = chrom2
    item.pos2 = pos2
    item.node_name = node_name
    item.read_enum = read_enum
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
                                      ReadEnum_t read_enum) nogil:

        cdef int idx, i, count_back, steps, node_name2
        cdef int sep = 0
        cdef int sep2 = 0
        cdef vector[int] found2
        cdef vector[int] found_exact
        cdef cpp_map[int, LocalVal]* forward_scope = &self.chrom_scope[chrom2]
        cdef cpp_map[int, LocalVal].iterator itr
        cdef cpp_pair[int, LocalVal] vitem

        # Re-initialize empty
        if current_chrom != self.local_chrom:
            self.local_chrom = current_chrom
            self.empty_scopes()

        if not self.loci.empty():  # Erase items out of range in local scope

            local_it = self.loci.lower_bound(current_pos - self.clst_dist)
            if local_it != self.loci.begin():
                self.loci.erase(self.loci.begin(), local_it)

            # Debug
            # local_it = forward_scope.begin()
            # while local_it != forward_scope.end():
            #     vitem = dereference(local_it)
            #     echo(vitem.first, vitem.second)
            #     preincrement(local_it)

            # Search FORWARD scope
            local_it = forward_scope.lower_bound(pos2)
            steps = 0
            if local_it != forward_scope.end():
                while steps < 4:
                    vitem = dereference(local_it)

                    if (read_enum == DELETION and vitem.second.read_enum == INSERTION) or (read_enum == INSERTION and vitem.second.read_enum == DELETION):
                    # Below: this increases speceficity slightly but hurts sensitivity for short-reads (del), but decreases both for insertions
                    # if (read_enum == DELETION and vitem.second.read_enum != DELETION) or (read_enum == INSERTION and vitem.second.read_enum != INSERTION) \
                    #       or (vitem.second.read_enum == DELETION and read_enum != DELETION) or (vitem.second.read_enum == INSERTION and read_enum != INSERTION):
                        steps += 1
                        continue  # Dont connect insertions to deletions

                    node_name2 = vitem.second.node_name
                    if node_name2 != node_name:  # Can happen due to within-read events
                        # if node_name == 9:
                        # echo("-->", node_name, node_name2, current_chrom == chrom2, current_pos, pos2, pos2 - current_pos, vitem.second.pos2 - vitem.first,
                        #          is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2),
                        #          )
                        #     echo(read_enum, vitem.second.read_enum, read_enum == INSERTION, vitem.second.read_enum == DELETION)
                        if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):

                            sep = c_abs(vitem.first - pos2)
                            sep2 = c_abs(vitem.second.pos2 - current_pos)
                            # echo("mad dist", self.max_dist, sep, sep2)
                            if sep < self.max_dist and sep2 < self.max_dist:
                                if sep < 35:
                                    found_exact.push_back(node_name2)
                                elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end):
                                    found2.push_back(node_name2)
                            elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end):
                                found2.push_back(node_name2)

                        elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end):
                            found2.push_back(node_name2)

                        # echo("")
                    if sep >= self.max_dist:
                        break  # No more to find

                    preincrement(local_it)
                    steps += 1
                    if local_it == forward_scope.end():
                        break

            if found_exact.empty():
                # Search lower
                local_it = forward_scope.lower_bound(pos2)
                if local_it != forward_scope.begin():
                    predecrement(local_it)  # Move back one before staring search, otherwise same value is processed twice

                    steps = 0
                    while steps < 4:
                        vitem = dereference(local_it)
                        if (read_enum == DELETION and vitem.second.read_enum == INSERTION) or (read_enum == INSERTION and vitem.second.read_enum == DELETION):
                        # if (read_enum == DELETION and vitem.second.read_enum != DELETION) or (read_enum == INSERTION and vitem.second.read_enum != INSERTION) \
                        #   or (vitem.second.read_enum == DELETION and read_enum != DELETION) or (vitem.second.read_enum == INSERTION and read_enum != INSERTION):
                            steps += 1
                            continue  # Dont connect insertions to deletions
                        node_name2 = vitem.second.node_name

                        if node_name2 != node_name:
                            # echo("2", node_name, node_name2, read_enum, vitem.second.clip_or_wr, is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2), span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2))
                            # if node_name == 10:
                            #     echo(current_pos, pos2, vitem.first, vitem.second.pos2)
                            if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):
                                sep = c_abs(vitem.first - pos2)
                                sep2 = c_abs(vitem.second.pos2 - current_pos)

                                if sep < self.max_dist and vitem.second.chrom2 == chrom2 and \
                                        sep2 < self.max_dist:
                                    if sep < 35:
                                        found_exact.push_back(node_name2)
                                    elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end):
                                        found2.push_back(node_name2)

                                elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end):
                                    found2.push_back(node_name2)

                            elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2, self.norm, self.thresh, read_enum, self.paired_end):
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
                       ReadEnum_t read_enum) nogil:

        # Add to scopes, if event is within read, add two references to forward scope. Otherwise when the
        # first break point drops from the local scope, the forward scope break will not connect with other
        # events that are added later. This is a fix for long read deletions mainly
        cdef cpp_map[int, LocalVal]* forward_scope  # = &self.chrom_scope[chrom2]
        if chrom2 == 10000000:
            forward_scope = &self.chrom_scope.back()  #[self.end_index]
        else:
            forward_scope = &self.chrom_scope[chrom2]

        # Add to local scope
        local_it = self.loci.find(current_pos)
        if local_it == self.loci.end():
            self.loci[current_pos] = make_local_val(chrom2, pos2, node_name, read_enum)
        else:
            dereference(local_it).second = make_local_val(chrom2, pos2, node_name, read_enum)

        # Add to forward scope
        local_it = forward_scope.find(pos2)
        if local_it == forward_scope.end():
            dereference(forward_scope)[pos2] = make_local_val(current_chrom, current_pos, node_name, read_enum)
        else:
            dereference(local_it).second = make_local_val(current_chrom, current_pos, node_name, read_enum)

        if read_enum == DELETION:
            local_it = forward_scope.find(current_pos)
            if local_it == forward_scope.end():
                dereference(forward_scope)[current_pos] = make_local_val(chrom2, pos2, node_name, read_enum)
            else:
                dereference(local_it).second = make_local_val(chrom2, pos2, node_name, read_enum)


cdef class TemplateEdges:
    cdef unordered_map[string, vector[int]] templates_s  # Better memory efficiency than dict, robin map was buggy for itertaing
    def __init__(self):
        pass

    cdef void add(self, str template_name, int flag, int node, int query_start):
        cdef vector[int] val
        cdef bytes key = bytes(template_name, encoding="utf8")
        # More efficient way of doing this?
        val.push_back(query_start)
        val.push_back(node)
        val.push_back(flag)
        self.templates_s[key].insert(self.templates_s[key].end(), val.begin(), val.end())

    def iterate_map(self):

        cdef unordered_map[string, vector[int]].iterator it = self.templates_s.begin()
        cdef string first
        cdef vector[int] second
        while it != self.templates_s.end():
            first = dereference(it).first
            second = dereference(it).second
            yield str(dereference(it).first), list(dereference(it).second)  # Array values are flag, node name, query start
            postincrement(it)


cdef void add_template_edges(G, TemplateEdges template_edges):
    # this function joins up template reads (read 1, read 2, plus any supplementary)
    cdef int ii, u_start, v_start, u, v, uflag, vflag
    # normally 2 reads for paired end, or >2 if supplementary reads
    for qname, arr in template_edges.iterate_map():
        read1_aligns = []
        read2_aligns = []
        for ii in range(0, len(arr), 3):
            if arr[ii + 2] & 64:
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
                    if not uflag & 2304:  # Is primary
                        primary1 = u
                    v_start, v, vflag = read1_aligns[ii + 1]
                    if not G.hasEdge(u, v):
                        G.addEdge(u, v, w=1)
                if primary1 is None:  # Check last in list
                    if not read1_aligns[-1][2] & 2304:
                        primary1 = read1_aligns[-1][1]

        if len(read2_aligns) > 0:
            if len(read2_aligns) == 1:
                if not read2_aligns[0][2] & 2304:  # Is primary
                    primary2 = read2_aligns[0][1]
            else:
                if len(read2_aligns) > 2:
                    read2_aligns = sorted(read2_aligns)
                for ii in range(len(read2_aligns) - 1):
                    u_start, u, uflag = read2_aligns[ii]
                    if not uflag & 2304:  # Is primary
                        primary2 = u
                    v_start, v, vflag = read2_aligns[ii + 1]
                    if not G.hasEdge(u, v):
                        G.addEdge(u, v, w=1)
                if primary2 is None:  # Check last in list
                    if not read2_aligns[-1][2] & 2304:
                        primary2 = read2_aligns[-1][1]

        if primary1 is not None and primary2 is not None:
            if not G.hasEdge(primary1, primary2):
                G.addEdge(primary1, primary2, w=1)


cdef class NodeToName:
    # Index these vectors to get the unique 'template_name'
    cdef vector[uint64_t] h
    cdef vector[uint16_t] f
    cdef vector[uint32_t] p
    cdef vector[uint16_t] c
    cdef vector[uint64_t] t
    cdef vector[int32_t] cigar_index
    cdef vector[uint32_t] event_pos

    def __cinit__(self):  # Possibly use a vector of structs instead. Reason for doing this way was to save mem
        # node names have the form (hash qname, flag, pos, chrom, tell, cigar index, event pos)
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
        return self.h[idx], self.f[idx], self.p[idx], self.c[idx], self.t[idx], self.cigar_index[idx], self.event_pos[idx]


cpdef tuple get_query_pos_from_cigarstring(cigar, pos):
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
        elif opp == "M":
            end += slen
            ref_end += slen
        elif opp == "D":
            ref_end += slen
        elif opp == "I":
            end += slen
        i = 1
    return start, end, pos, ref_end


cpdef tuple get_query_pos_from_cigartuples(r):
    # Infer the position on the query sequence of the alignment using cigar string
    cdef int end = 0
    cdef int start = 0
    cdef bint i = 0
    cdef int ref_end = r.pos
    cdef int opp, slen
    for opp, slen in r.cigartuples:
        if opp == 0:
            ref_end += slen
            end += slen
        elif opp == 4 or opp == 5:
            if not i:
                start += slen
                end += slen
                i = 1
            else:
                break
        elif opp == 1:
            end += slen
        elif opp == 2:  # deletion
            ref_end += slen
        i = 1
    return start, end, r.pos, ref_end, r.rname, r.mapq, 0


cdef alignments_from_sa_tag(r, gettid, thresh, paired_end, mapq_thresh):
    # Puts other alignments in order of query position, gets index of the current alignment
    cdef int query_length = r.infer_read_length()  # Note, this also counts hard-clips
    cdef int chrom2, start_pos2, query_start, query_end, ref_start, ref_end, start_temp, aln_start, aln_end

    current_strand = "-" if r.flag & 16 else "+"
    query_aligns = [get_query_pos_from_cigartuples(r)]
    aln_start = query_aligns[0][2]
    aln_end = query_aligns[0][3]
    aln_chrom = r.rname

    for sa_block in r.get_tag("SA").split(";"):
        if sa_block == "":
            break  # End

        sa = sa_block.split(",", 5)
        mq = int(sa[4])  # filter by mapq at later stage
        chrom2 = gettid(sa[0])
        start_pos2 = int(sa[1])
        strand = sa[2]
        cigar = sa[3]
        matches = [(int(slen), opp) for slen, opp in re.findall(r'(\d+)([A-Z]{1})', sa[3])]  # parse cigar
        query_start, query_end, ref_start, ref_end = get_query_pos_from_cigarstring(matches, start_pos2)

        if current_strand != strand:  # count from end
            start_temp = query_length - query_end
            query_end = start_temp + query_end - query_start
            query_start = start_temp

        # If another local alignment is found use only this, usually corresponds to the other side of an insertion/dup
        if aln_chrom == chrom2 and position_distance(aln_start, aln_end, ref_start, ref_end) < thresh:
            query_aligns = [query_aligns[0], (query_start, query_end, ref_start, ref_end, chrom2, mq)]
            break

        query_aligns.append((query_start, query_end, ref_start, ref_end, chrom2, mq))

    query_aligns = sorted(query_aligns)
    cdef int index = 0
    for index in range(len(query_aligns)):
        if len(query_aligns[index]) == 7:
            break

    return query_aligns, index


cdef connect_right(a, b, r, paired, max_dist, mq_thresh):
    event_pos = a[3]
    chrom = a[4]
    pos2 = b[2]
    chrom2 = b[4]
    pos2_mq = b[5]
    if pos2_mq < mq_thresh:
        return event_pos, chrom, event_pos, chrom, ReadEnum_t.BREAKEND
    if paired:
        # If paired, and SA block fits between two normal primary alignments, and block is not local, then
        # ignore block and try and call insertion
        if not r.flag & 2304 and \
            r.flag & 2 and r.pnext >= event_pos and \
            (chrom != chrom2 or abs(pos2 - event_pos) > max_dist):  # is primary, is proper pair. Use primary mate info, not SA
            return event_pos, chrom, event_pos, chrom, ReadEnum_t.BREAKEND

    return event_pos, chrom, pos2, chrom2, 0


cdef connect_left(a, b, r, paired, max_dist, mq_thresh):
    event_pos = a[2]
    chrom = a[4]
    pos2 = b[3]
    chrom2 = b[4]
    pos2_mq = b[5]
    if pos2_mq < mq_thresh:
        return event_pos, chrom, event_pos, chrom, ReadEnum_t.BREAKEND
    if paired:
        if not r.flag & 2304 and \
            r.flag & 2 and r.pnext <= event_pos and \
            (chrom != chrom2 or abs(pos2 - event_pos) > max_dist):
            return event_pos, chrom, event_pos, chrom, ReadEnum_t.BREAKEND

    return event_pos, chrom, pos2, chrom2, 0


cdef cluster_clipped(G, r, ClipScoper_t clip_scope, chrom, pos, node_name):

    cdef int other_node
    cdef unordered_set[int] clustered_nodes
    clip_scope.update(r, node_name, chrom, pos, clustered_nodes)
    if not clustered_nodes.empty():
        for other_node in clustered_nodes:

            if not G.hasEdge(node_name, other_node):
                G.addEdge(node_name, other_node, 2)


cdef bint add_to_graph(G, AlignedSegment r, PairedEndScoper_t pe_scope, TemplateEdges_t template_edges,
                       NodeToName node_to_name, genome_scanner,
                       int flag, int chrom, tell, int cigar_index, int event_pos,
                       int chrom2, int pos2, ClipScoper_t clip_scope, ReadEnum_t read_enum,
                       bint p1_overlaps, bint p2_overlaps, bint mm_only):

    # Adds relevant information to graph and other data structures for further processing
    cdef vector[int] other_nodes  # Other alignments to add edges between
    cdef int node_name = G.addNode()
    cdef uint64_t v = xxhasher(bam_get_qname(r._delegate), len(r.qname), 42)  # Hash qname to save mem

    node_to_name.append(v, flag, r.pos, chrom, tell, cigar_index, event_pos)  # Index this list to get the template_name
    genome_scanner.add_to_buffer(r, node_name, tell)  # Add read to buffer

    if read_enum < 2:  # Prevents joining up within-read svs with between-read svs
        template_edges.add(r.qname, flag, node_name, r.query_alignment_start)

    both_overlap = p1_overlaps and p2_overlaps

    if read_enum != BREAKEND and (not mm_only or not both_overlap):
        other_nodes = pe_scope.find_other_nodes(node_name, chrom, event_pos, chrom2, pos2, read_enum)

    elif chrom != chrom2:  # Note all BREAKENDS have chrom != chrom2, but also includes translocations or read_enum == BREAKEND:
        cluster_clipped(G, r, clip_scope, chrom, event_pos, node_name)

        # if read_enum == BREAKEND:  # Try and connect soft-clips to within-read events (doing this hurts precision)
        #     other_nodes = pe_scope.find_other_nodes(node_name, chrom, event_pos, chrom, event_pos, read_enum)
    if not mm_only or not both_overlap:
        pe_scope.add_item(node_name, chrom, event_pos, chrom2, pos2, read_enum)
    # Debug:
    # echo("---", r.qname, read_enum, node_name, event_pos, pos2, list(other_nodes), chrom, chrom2)
    # look = set(['V300082976L4C001R0311226430'])
    # node_look = set([659281, 659282, 659283, 659284, 659285, 659286])
    # if r.qname in look:
    # if r.qname == "HISEQ2500-10:539:CAV68ANXX:7:2214:16836:93172":
    #     echo("@", r.flag, node_name, event_pos, pos2, list(other_nodes), cigar_index)
    #     echo(both_overlap, "enum", read_enum, p1_overlaps, p2_overlaps)
    #     quit()
    if not other_nodes.empty():
        for other_node in other_nodes:
            if not G.hasEdge(node_name, other_node):
                # 2 signifies that this is a local edge as opposed to a template edge
                G.addEdge(node_name, other_node, 2)
        return 1
    return 0


cdef int good_quality_clip(AlignedSegment r, int clip_length):

    # Use sliding window to check that soft-clip has at least n good bases over qual 20
    cdef const unsigned char[:] quals = r.query_qualities
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
                        # Make sure stretch is not homopolymer
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


cdef void process_alignment(G, AlignedSegment r, int clip_l, int loci_dist, gettid,
                            overlap_regions, int clustering_dist, PairedEndScoper_t pe_scope,
                            int cigar_index, int event_pos, int paired_end, long tell, genome_scanner,
                            TemplateEdges_t template_edges, NodeToName node_to_name,
                            int cigar_pos2, int mapq_thresh, ClipScoper_t clip_scope,
                            ReadEnum_t read_enum, bad_clip_counter, bint mm_only):

    # Determines where the break point on the alignment is before adding to the graph
    cdef int other_node, clip_left, clip_right
    cdef int current_overlaps_roi, next_overlaps_roi
    cdef bint add_primark_link, add_insertion_link
    cdef int chrom = r.rname
    cdef int chrom2 = r.rnext
    cdef int pos2 = r.pnext
    cdef int flag = r.flag
    cdef int left_clip, right_clip
    cdef str qname = r.qname
    cdef ReadEnum_t inferred_clip_type
    cdef uint64_t v

    if paired_end and read_enum == SPLIT and flag & 8:  # clip event, or whole read, but mate is unmapped
        return

    cdef bint success
    cdef bint good_clip

    if paired_end or flag & 1:

        # Cluster paired-end mates
        pnext = r.pnext
        rnext = r.rnext
        chrom2 = rnext
        pos2 = pnext

        current_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, r.rname, event_pos, event_pos+1)
        next_overlaps_roi = False
        if current_overlaps_roi:  # check if both sides of SV are in overlaps regions
            next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, chrom2, pos2, pnext+1)

        # Special treatment of supplementary and local reads; need to decide where the partner is
        # Either use the rnext:pnext or a site(s) listed in the SA tag: The rnext:pext can often be at the
        # same loci as the read which leads to problems when linking w=2 edges

        add_primary_link = 1
        add_insertion_link = 0

        good_clip = good_quality_clip(r, 15)
        if not good_clip:
            bad_clip_counter.add(chrom, r.pos)

        if read_enum == SPLIT:
            # Parse SA tag. For paired reads
            if r.has_tag("SA") and good_clip:  # Parse SA, first alignment is the other read primary line
                all_aligns, index = alignments_from_sa_tag(r, gettid, loci_dist, paired_end, mapq_thresh)
                event = all_aligns[index]
                if len(all_aligns) == 1:
                    return  # should'nt happen

                if index < len(all_aligns) - 1:  # connect to next
                    event_pos, chrom, pos2, chrom2, inferred_clip_type = connect_right(all_aligns[index], all_aligns[index + 1], r, paired_end, loci_dist, mapq_thresh)
                    cigar_index = len(r.cigartuples) - 1
                    if inferred_clip_type == BREAKEND:
                        read_enum = ReadEnum_t.BREAKEND

                if index > 0:
                    event_pos, chrom, pos2, chrom2, inferred_clip_type = connect_left(all_aligns[index], all_aligns[index -1], r, paired_end, loci_dist, mapq_thresh)
                    cigar_index = 0
                    if inferred_clip_type == BREAKEND:
                        read_enum = ReadEnum_t.BREAKEND

        elif read_enum == DISCORDANT:
            pass  # rnext and pnext set as above

        if read_enum == BREAKEND:
            if mm_only and current_overlaps_roi and next_overlaps_roi:
                # Probably too many reads in ROI to reliably separate out break end reads
                return
            if good_clip and cigar_clip(r, clip_l):
                read_enum = BREAKEND
                chrom2 = 10000000
                pos2 = event_pos

        if read_enum == DELETION or read_enum == INSERTION:  # DELETION or INSERTION within read
            chrom2 = chrom
            if r.cigartuples[cigar_index][0] != 1:  # not insertion, use length of cigar event
                pos2 = cigar_pos2
            else:
                pos2 = event_pos

        success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                               tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum, current_overlaps_roi, next_overlaps_roi,
                               mm_only)

    ###
    else:  # Single end
        current_overlaps_roi, next_overlaps_roi = False, False  # not supported yet
        if read_enum == SPLIT:
            # Use SA tag to get chrom2 and pos2
            if r.has_tag("SA"):

                all_aligns, index = alignments_from_sa_tag(r, gettid, loci_dist, paired_end, mapq_thresh)
                event = all_aligns[index]

                if len(all_aligns) == 1:
                    return  # shouldnt happen

                if index < len(all_aligns) - 1:  # connect to next
                    event_pos, chrom, pos2, chrom2, _ = connect_right(all_aligns[index], all_aligns[index + 1], r, paired_end, loci_dist, mapq_thresh)
                    cigar_index = 0
                    success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                           tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum,
                                           current_overlaps_roi, next_overlaps_roi,
                                           mm_only)

                if index > 0:
                    event_pos, chrom, pos2, chrom2, _ = connect_left(all_aligns[index], all_aligns[index -1], r, paired_end, loci_dist, mapq_thresh)
                    cigar_index = len(r.cigartuples) - 1
                    success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                           tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum,
                                           current_overlaps_roi, next_overlaps_roi,
                                           mm_only)
        elif read_enum >= 2:  # Sv within read
            chrom2 = r.rname
            if r.cigartuples[cigar_index][0] != 1:  # If not insertion
                pos2 = cigar_pos2  # event_pos + r.cigartuples[cigar_index][1]
            else:
                pos2 = event_pos
            success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                   tell, cigar_index, event_pos, chrom2, pos2, clip_scope, read_enum,
                                   current_overlaps_roi, next_overlaps_roi, mm_only)


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


cpdef tuple construct_graph(genome_scanner, infile, int max_dist, int clustering_dist, int k=16, int m=7, int clip_l=21,
                            int min_sv_size=30,
                            int minimizer_support_thresh=2, int minimizer_breadth=3,
                            int minimizer_dist=10, int mapq_thresh=1, debug=None, procs=1,
                            int paired_end=1, int read_length=150, bint contigs=True,
                            float norm_thresh=100, float spd_thresh=0.3, bint mm_only=False):

    t0 = time.time()
    logging.info("Building graph with clustering distance {} bp, scope length {} bp".format(max_dist, clustering_dist))

    cdef TemplateEdges_t template_edges = TemplateEdges()  # Edges are added between alignments from same template, after building main graph
    cdef int event_pos, cigar_index, opp, length
    node_to_name = NodeToName()  # Map of nodes -> read ids


    cdef ClipScoper_t clip_scope = ClipScoper(minimizer_dist, k=k, m=m, clip_length=clip_l,  # Keeps track of local reads
                       minimizer_support_thresh=minimizer_support_thresh,
                       minimizer_breadth=minimizer_breadth, read_length=read_length)


    # Infers long-range connections, outside local scope using pe information
    cdef PairedEndScoper_t pe_scope = PairedEndScoper(max_dist, clustering_dist, infile.header.nreferences, norm_thresh, spd_thresh, paired_end)

    bad_clip_counter = BadClipCounter(infile.header.nreferences)

    cdef long tell
    G = map_set_utils.Py_SimpleGraph()
    overlap_regions = genome_scanner.overlap_regions  # Get overlapper, intersect reads with intervals
    gettid = infile.gettid

    cdef int pos2
    cdef bint added
    cdef ReadEnum_t read_enum

    # ReadEnum_t read_enum
    # an indicator to describe the kind of sv/alignment
    # 1 = soft-clip or a split-read, the supplementary mapping might be discarded. whole read is a node
    # 2 = deletion event contained within read
    # 3 = insertion event contained within read
    # 4 = mate-pair is unmapped, or supplementary uninformative, probably an insertion

    cdef bint clipped
    cdef vector[CigarEvent] events_to_add
    cdef vector[CigarEvent].iterator itr_events
    cdef CigarEvent v

    for chunk in genome_scanner.iter_genome():
        for r, tell in chunk:
            if r.mapq < mapq_thresh:
                continue

            pos2 = -1
            event_pos = r.pos
            added = False
            clipped = 0
            events_to_add.clear()

            if len(r.cigartuples) > 1:

                if r.has_tag("SA"):

                    # Set cigar-index to -1 means it is unset, will be determined during SA parsing
                    cigar_index = -1
                    process_alignment(G, r, clip_l, max_dist, gettid,
                                      overlap_regions, clustering_dist, pe_scope,
                                      cigar_index, event_pos, paired_end, tell, genome_scanner,
                                      template_edges, node_to_name, pos2, mapq_thresh, clip_scope, ReadEnum_t.SPLIT,
                                      bad_clip_counter, mm_only)
                    added = True

                for cigar_index, (opp, length) in enumerate(r.cigartuples):

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
                        elif opp == 4:
                            clipped = 1

            if not added and contigs:
                # Whole alignment will be used, try infer position from soft-clip
                cigar_index = -1
                pos2 = -1
                left_clip_size, right_clip_size = clip_sizes_hard(r)  # soft and hard-clips
                if r.flag & 8 and clipped:
                    # skip if both ends are clipped, usually means its a chunk of badly mapped sequence
                    if not (left_clip_size and right_clip_size) and good_quality_clip(r, 20):
                        # Mate is unmapped, insertion type. Only add if soft-clip is available
                        process_alignment(G, r, clip_l, max_dist, gettid,
                                          overlap_regions, clustering_dist, pe_scope,
                                          cigar_index, event_pos, paired_end, tell, genome_scanner,
                                          template_edges, node_to_name,
                                          pos2, mapq_thresh, clip_scope, ReadEnum_t.BREAKEND, bad_clip_counter,
                                          mm_only)
                else:
                    # Use whole read, could be normal or discordant
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
                                      mm_only)

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
                                      mm_only)

                    preincrement(itr_events)

    add_template_edges(G, template_edges)

    return G, node_to_name, bad_clip_counter


cpdef dict get_reads(infile, sub_graph_reads):

    rd = dict()
    cdef int j, int_node
    cdef long int p
    cdef uint64_t v
    cdef AlignedSegment a
    for int_node, node in sub_graph_reads.items():
        node = tuple(node[:-2])  # drop cigar index and event pos
        p = node[4]
        infile.seek(p)
        a = next(infile)
        v = xxhasher(bam_get_qname(a._delegate), len(a.qname), 42)
        n1 = (v, a.flag, a.pos, a.rname, p)
        # Try next few reads, sometimes they are on top of one another
        if n1 != node:
            for j in range(5):
                a = next(infile)
                n2 = (xxhasher(bam_get_qname(a._delegate), len(a.qname), 42), a.flag, a.pos, a.rname, p)
                if n2 == node:
                    rd[int_node] = a
                    break
        else:
            rd[int_node] = a
    return rd


cdef set BFS_local(G, int source, unordered_set[int]& visited ):
    # Create a queue for BFS
    queue = array.array("L", [source])  # [source]  # consider vector
    nodes_found = set([])
    cdef int u, v
    while queue:
        u = queue.pop(0)
        for v in G.neighbors(u):
            if visited.find(v) == visited.end(): #v not in visited:
                if G.weight(u, v) > 1:
                    if u not in nodes_found:
                        nodes_found.add(u)
                    if v not in nodes_found:
                        nodes_found.add(v)
                        queue.append(v)
        visited.insert(u)

    return nodes_found


cdef dict get_partitions(G, nodes):

    cdef unordered_set[int] seen
    cdef int u, v, i
    parts = []
    for u in nodes:
        if seen.find(u) != seen.end(): #u in seen:
            continue

        for v in G.neighbors(u):
            if seen.find(v) != seen.end(): #v in seen:
                continue

            if G.weight(u, v) > 1:  # weight is 2 or 3, for normal or black edges
                found = BFS_local(G, u, seen)
                if found:
                    parts.append(array.array("L", found))

        seen.insert(u)
    return {i: p for i, p in enumerate(parts)}


cdef tuple count_support_between(G, parts, int min_support):

    cdef int i, j, node, child, any_out_edges
    cdef tuple t

    if len(parts) == 0:
        return {}, {}
    elif len(parts) == 1:
        return {}, {list(parts.keys())[0]: array.array("L", list(parts.values())[0])}

    # Make a table to count from, int-int
    cdef Py_Int2IntMap p2i = map_set_utils.Py_Int2IntMap()
    for i, p in parts.items():
        for node in p:
            p2i.insert(node, i)

    # Count the links between partitions. Split reads into sets ready for calling
    # No counting of read-pairs templates or 'support', just a count of linking alignments
    # counts (part_a, part_b): {part_a: {node 1, node 2 ..}, part_b: {node4, ..} }
    counts = {}
    self_counts = {}

    seen_t = set([])
    for i, p in parts.items():

        current_t = set([])
        for node in p:
            any_out_edges = 0  # Keeps track of number of outgoing pairs, or self edges

            for child in G.neighbors(node):

                if not p2i.has_key(child):
                    continue  # Exterior child, not in any partition

                # Partition of neighbor node
                j = p2i.get(child)
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

    return counts, self_counts


cpdef dict proc_component(node_to_name, component, read_buffer, infile, G,
                         int min_support, int procs, int paired_end):

    n2n = {}
    reads = {}
    cdef int support_estimate = 0
    cdef int v
    for v in component:

        if procs == 1 and v in read_buffer:
            reads[v] = read_buffer[v]
        # Need to keep a record of all node info, and cigar indexes
        key = node_to_name[v]
        if key[5] != -1:
            support_estimate += 2
        else:
            support_estimate += 1
        n2n[v] = key

    if support_estimate < min_support:
        return {}

    # Explore component for locally interacting nodes; create partitions using these
    partitions = get_partitions(G, component)
    support_between, support_within = count_support_between(G, partitions, min_support)

    if len(support_between) == 0 and len(support_within) == 0:
        if not paired_end:

            if len(n2n) >= min_support or len(reads) >= min_support:
                return {"parts": {}, "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
            else:
                return {}

        else:
            # single paired end template can have 3 nodes e.g. two reads plus supplementary
            if min_support == 1 and (len(n2n) >= min_support or len(reads) >= min_support):
                return {"parts": {}, "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
            elif len(reads) >= min_support:
                return {"parts": {}, "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
            else:
                return {}

    sb = {}
    for edge, vd in support_between.items():
        sb[edge] = vd

    # Debug:
    # echo("parts", partitions)
    # echo("s_between", sb)
    # echo("s_within", support_within)
    # quit()
    # echo("n2n", n2n.keys())
    # node_look = set(range(653526, 653532))
    # node_look = set(range(8))
    # echo(node_look.intersection(set(n2n.keys())))
    # if node_look.intersection(set(n2n.keys())):
    #     echo("parts", partitions)
    #     echo("s_between", sb)
    #     echo("s_within", support_within)

    return {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}
