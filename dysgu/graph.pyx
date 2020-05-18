#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

from __future__ import absolute_import
import time
import click
from collections import defaultdict, deque

# import mmh3
import ncls

import numpy as np
cimport numpy as np
from cpython cimport array
import array
import re

from libcpp.vector cimport vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string

from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t
from libc.stdlib cimport abs as c_abs

from cython.operator import dereference, postincrement, postdecrement, preincrement, predecrement

from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport bam_get_qname

from dysgu cimport map_set_utils
from dysgu import io_funcs

from dysgu.map_set_utils cimport unordered_set
from dysgu.map_set_utils cimport hash as xxhasher

ctypedef cpp_pair[int, int] cpp_item

ctypedef map_set_utils.Py_IntSet Py_IntSet
ctypedef map_set_utils.Py_Int2IntMap Py_Int2IntMap

ctypedef cpp_map[int, cpp_item] ScopeItem_t
ctypedef vector[ScopeItem_t] ForwardScope_t

ctypedef cpp_pair[int, cpp_item] event_item
ctypedef long int long_int
ctypedef PairedEndScoper PairedEndScoper_t
ctypedef TemplateEdges TemplateEdges_t

ctypedef NodeToName NodeToName_t

def echo(*args):
    click.echo(args, err=True)


cdef class Table:
    cdef vector[np.int64_t] starts
    cdef vector[np.int64_t] ends
    cdef vector[np.int64_t] values

    cpdef void add(self, int s, int e, int v):
        # with nogil:
            self.starts.push_back(s)
            self.ends.push_back(e)
            self.values.push_back(v)

    def get_val(self, v):
        cdef vector[np.int64_t] values = v
        cdef np.ndarray[np.int64_t] a = np.empty(values.size(), dtype=np.int)
        cdef int len_a = len(a)  #.shape[0]
        cdef int i
        with nogil:
            for i in range(len_a):
                a[i] = values[i]
        return a

    def containment_list(self):
        return ncls.NCLS(self.get_val(self.starts), self.get_val(self.ends), self.get_val(self.values))


cdef set sliding_window_minimum(int k, int m, str s):
    """End minimizer. A iterator which takes the size of the window, `k`, and an iterable,
    `li`. Then returns an iterator such that the ith element yielded is equal
    to min(list(li)[max(i - k + 1, 0):i+1]).
    Each yield takes amortized O(1) time, and overall the generator takes O(k)
    space.
    https://github.com/keegancsmith/Sliding-Window-Minimum/blob/master/sliding_window_minimum.py"""

    # Cpp version
    cdef int i = 0
    cdef int end = len(s) - m + 1
    cdef cpp_deque[cpp_item] window2
    cdef long int hx2

    cdef bytes s_bytes = s.encode("ascii")
    # cdef cpp_set[int] seen2  # Using
    # cdef cpp_u_set[int] seen2
    seen2 = set([])

    # cdef cpp_item last
    for i in range(end):
        # xxhasher(bam_get_qname(r._delegate), len(qname), 42)
        # hx2 = mmh3.hash(s[i:i+m], 42)
        # s_bytes = s[i:i+m].encode("ascii")
        hx2 = xxhasher(s_bytes[i:i+m], len(s_bytes), 42)
        while window2.size() != 0 and window2.back().first >= hx2:
            window2.pop_back()

        window2.push_back(cpp_item(hx2, i))
        while window2.front().second <= i - k:
            window2.pop_front()

        i += 1

        minimizer_i = window2.front().first

        #if minimizer_i not in seen2:
        seen2.add(minimizer_i)
        # if seen2.find(minimizer_i) == seen2.end():
        #     seen2.insert(minimizer_i)

    return seen2  #set(seen2)


cdef str left_soft_clips(str seq, int code_length):
    return seq[0:code_length]


cdef str right_soft_clips(str seq, int code_length):
    return seq[len(seq) - code_length:]


class ClipScoper:
    """Keeps track of which reads are in scope. Maximum distance depends on the template insert_median"""
    def __init__(self, int max_dist, int k, int m, int clip_length, int minimizer_support_thresh,
                 int minimizer_breadth):
        self.max_dist = max_dist
        self.k = k
        self.m = m
        self.clip_length = clip_length
        self.scope = deque([])
        self.minimizer_table = defaultdict(set) # For left and right clips, reference counter
        self.read_minimizers = defaultdict(set)
        self.current_chrom = 0
        self.minimizer_support_thresh = minimizer_support_thresh
        self.minimizer_matches = minimizer_breadth

    def _add_m_find_candidates(self, clip_set, int name, int idx):
        # idx 0 is for left clips, 1 is for right clips
        target_counts = defaultdict(int)
        cdef int find_candidate = 1
        cdef int m
        cdef int item
        res = set([])
        for m in clip_set:

            min_id = (m, idx)
            self.read_minimizers[name].add(min_id)

            if min_id not in self.minimizer_table:
                self.minimizer_table[min_id].add(name)
                continue

            elif find_candidate:  # Look for suitable partners
                targets = self.minimizer_table[min_id]

                for item in targets:

                    target_counts[item] += 1

                    if target_counts[item] >= self.minimizer_matches and \
                            len(target_counts) >= self.minimizer_support_thresh:
                        res.add(item)
                    if len(res) >= 3:  # Maximum of 3 edges for each read
                        find_candidate = 0
                        break

            self.minimizer_table[min_id].add(name)

        # Return best 3 edges, stop graph getting too dense
        return res

    def _insert(self, str seq, int cigar_start, int cigar_end, int input_read):
        # Find soft-clips of interest
        cdef set clip_set
        cdef str clip_seq
        targets = set([])

        if cigar_start >= self.clip_length:
            clip_seq = left_soft_clips(seq, cigar_start)
            clip_set = sliding_window_minimum(self.k, self.m, clip_seq)
            targets |= self._add_m_find_candidates(clip_set, input_read, 0)

        if cigar_end >= self.clip_length:
            clip_seq = right_soft_clips(seq, cigar_end)
            clip_set = sliding_window_minimum(self.k, self.m, clip_seq)
            targets |= self._add_m_find_candidates(clip_set, input_read, 1)

        return targets

    def _pop_minimizers(self):
        cdef tuple name_pos = self.scope.popleft()
        cdef int name = name_pos[0]
        cdef set t
        try:  # Read with same name was already processed, very rare but can happen
            t = self.read_minimizers[name]
        except KeyError:
            return
        cdef tuple item
        for item in self.read_minimizers[name]:
            try:
                self.minimizer_table[item].remove(name)
                if len(self.minimizer_table[item]) == 0:
                    del self.minimizer_table[item]
            except KeyError:
                pass
        del self.read_minimizers[name]

    def update(self, int input_read, str seq, int cigar_start, int cigar_end, int chrom, int position):

        if len(self.scope) == 0:
            self.scope.append((input_read, position))
            self.current_chrom = chrom
            self._insert(seq, cigar_start, cigar_end, input_read)  # Add minimizers
            return

        elif chrom != self.current_chrom:
            self.scope = deque([(input_read, position)])  # Empty scope on new chromosome
            self.minimizer_table = defaultdict(set)
            self.read_minimizers = defaultdict(set)
            self.current_chrom = chrom
            self._insert(seq, cigar_start, cigar_end, input_read)
            return

        # Remove out of scope reads and minimizers
        while True:
            if len(self.scope) == 0:
                break
            if abs(self.scope[0][1] - position) > self.max_dist:
                self._pop_minimizers()
            else:
                break
        self.scope.append((input_read, position))

        cdef set d = self._insert(seq, cigar_start, cigar_end, input_read)

        if len(d) > 0:
            return d  # max(d, key=d.get)  # Best candidate with most overlapping minimizers


cdef struct LocalVal:
    int chrom2
    int pos2
    int node_name
    uint8_t clip_or_wr


cdef LocalVal make_local_val(int chrom2, int pos2, int node_name, uint8_t clip_or_wr) nogil:
    cdef LocalVal item
    item.chrom2 = chrom2
    item.pos2 = pos2
    item.node_name = node_name
    item.clip_or_wr = clip_or_wr

    return item


# cdef bint is_overlapping(int x1, int x2, int y1, int y2) nogil:
#     return max(x1, y1) <= min(x2, y2)


cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil:
    # Insertions have same x1/y1 position
    if x1 == x2:  # This might confuse insertions with deletions, but links soft-clipped reads with indels
        if x1 == y1 or x1 == y2:
            return True
    if y1 == y2:
        if y1 == x1 or y1 == x2:
            return True

    cdef int temp_v
    if x2 < x1:
        temp_v = x2
        x2 = x1
        x1 = temp_v
    if y2 < y1:
        temp_v = y2
        y2 = y1
        y1 = temp_v
    cdef float overlap = float(max(0, (min(x2, y2) - max(x1, y1))))

    if overlap == 0:
        return False
    if (overlap / float(c_abs(x2 - x1))) > 0.3 and (overlap / float(c_abs(y2 - y1))) > 0.3:
        return True


cdef class PairedEndScoper:

    cdef int clst_dist
    cdef int max_dist
    cdef int local_chrom
    cdef cpp_map[int, LocalVal] loci  # Track the local breaks and mapping locations
    cdef vector[cpp_map[int, LocalVal]] chrom_scope  # Track the mate-pair breaks and locations

    def __init__(self, max_dist, clst_dist, n_references):
        self.clst_dist = clst_dist
        self.max_dist = max_dist # * 1.5
        self.local_chrom = -1

        cdef cpp_map[int, LocalVal] scope
        for n in range(n_references):
            self.chrom_scope.push_back(scope)

    cdef void empty_scopes(self) nogil:
        for idx in range(self.chrom_scope.size()):
            if not self.chrom_scope[idx].empty():
                self.chrom_scope[idx].clear()
        self.loci.clear()

    cdef vector[int] update(self, int node_name, int current_chrom, int current_pos, int chrom2, int pos2,
                            int clip_or_wr) nogil:

        cdef int idx, i, count_back, steps, node_name2 #, lowest_pos
        cdef int sep = 0
        cdef int sep2 = 0
        cdef vector[int] found2
        cdef vector[int] found_exact
        cdef cpp_map[int, LocalVal]* forward_scope = &self.chrom_scope[chrom2]
        cdef cpp_map[int, LocalVal].iterator itr
        cdef cpp_pair[int, LocalVal] vitem

        # echo("current node", node_name, current_chrom, current_pos, chrom2, pos2, forward_scope.size(), clip_or_wr)
        # Re-initialize empty
        if current_chrom != self.local_chrom:
            self.local_chrom = current_chrom
            self.empty_scopes()

        if not self.loci.empty():

            # local_it = forward_scope.begin()
            # while local_it != forward_scope.end():
            #     vitem = dereference(local_it)
            #     echo("before", vitem.first, vitem.second)
            #     preincrement(local_it)

            # Erase items out of range in local scope
            # lowest_pos = current_pos
            # if current_chrom == chrom2 and pos2 < current_pos:  # Make sure only clear up to lowest position
            #     lowest_pos = pos2

            local_it = self.loci.lower_bound(current_pos - self.clst_dist)
            if local_it != self.loci.begin():
                # echo("dropping up to", dereference(local_it).first)
                self.loci.erase(self.loci.begin(), local_it)

            if current_chrom != chrom2:  # and forward scope
                local_it = dereference(forward_scope).lower_bound(current_pos - self.clst_dist)
                if local_it != dereference(forward_scope).begin():
                    dereference(forward_scope).erase(dereference(forward_scope).begin(), local_it)

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

                    node_name2 = vitem.second.node_name
                    if node_name2 != node_name:  # Can happen due to within-read events

                        if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):
                            # echo(node_name2, is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2))

                            sep = c_abs(vitem.first - pos2)
                            sep2 = c_abs(vitem.second.pos2 - current_pos)

                            if sep < self.max_dist and vitem.second.chrom2 == chrom2 and \
                                    sep2 < self.max_dist:

                                if sep < 25 and (clip_or_wr > 0 or vitem.second.clip_or_wr):
                                    found_exact.push_back(node_name2)
                                else:
                                    found2.push_back(node_name2)

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
                        node_name2 = vitem.second.node_name
                        if node_name2 != node_name:
                            if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2): # is_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):

                                # echo(node_name2, is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2))

                                sep = c_abs(vitem.first - pos2)
                                sep2 = c_abs(vitem.second.pos2 - current_pos)
                                # if node_name == 55:
                                #     echo(sep, sep2, current_pos, pos2, vitem.second)
                                # echo("2node2", node_name2, sep, sep < self.max_dist, vitem.second.chrom2 == chrom2,
                                #         c_abs(vitem.second.pos2 - current_pos), self.max_dist, )

                                if sep < self.max_dist and vitem.second.chrom2 == chrom2 and \
                                        sep2 < self.max_dist:
                                    if sep < 25 and (clip_or_wr > 0 or vitem.second.clip_or_wr):
                                        found_exact.push_back(node_name2)
                                    else:
                                        found2.push_back(node_name2)

                        if local_it == forward_scope.begin() or sep >= self.max_dist:
                            break
                        predecrement(local_it)
                        steps += 1

        # Add tto scope
        if current_pos == pos2 and current_chrom == chrom2:  # Update same position for insertion
            local_it = forward_scope.find(pos2)
            if local_it == forward_scope.end():
                dereference(forward_scope)[pos2] = make_local_val(current_chrom, current_pos, node_name, clip_or_wr)
            else:
                dereference(local_it).second = make_local_val(current_chrom, current_pos, node_name, clip_or_wr)
        else:
            # Add to scopes, if event is within read, add two references to forward scope. Otherwise when the
            # first break point drops from the local scope, the forward scope break will not connect with other
            # events that are added later. This is a fix for long reads

            # Add to local scope
            local_it = self.loci.find(current_pos)
            if local_it == self.loci.end():
                self.loci[current_pos] = make_local_val(chrom2, pos2, node_name, clip_or_wr)
            else:
                dereference(local_it).second = make_local_val(chrom2, pos2, node_name, clip_or_wr)

            # Add to forward scope
            local_it = forward_scope.find(pos2)
            if local_it == forward_scope.end():
                dereference(forward_scope)[pos2] = make_local_val(current_chrom, current_pos, node_name, clip_or_wr)
            else:
                dereference(local_it).second = make_local_val(current_chrom, current_pos, node_name, clip_or_wr)

            if clip_or_wr == 2:
                local_it = forward_scope.find(current_pos)
                if local_it == forward_scope.end():
                    dereference(forward_scope)[current_pos] = make_local_val(chrom2, pos2, node_name, clip_or_wr)
                else:
                    dereference(local_it).second = make_local_val(chrom2, pos2, node_name, clip_or_wr)

        # echo(found_exact, found2)
        # if node_name == 96:
        #     quit()
        if not found_exact.empty():
            return found_exact
        else:
            return found2


cdef class TemplateEdges:

    cdef unordered_map[string, vector[int]] templates_s  # Better memory efficiency than dict -> use robinmap?

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
            yield str(dereference(it).first), list(dereference(it).second)
            postincrement(it)
            # Array values are flag, node name, query start


cdef void add_template_edges(G, TemplateEdges template_edges): #TemplateEdges template_edges):

    cdef int ii, u_start, v_start, u, v, uflag, vflag

    # for qname, arr in template_edges.templates.items():  # normally 2 reads, or >2 if supplementary reads
    # for qname, (read1_aligns, read2_aligns) in template_edges.templates.items():
    for qname, arr in template_edges.iterate_map():
        # echo(qname, arr)
        read1_aligns = []
        read2_aligns = []
        for ii in range(0, len(arr), 3):
            if arr[ii + 2] & 64:
                read1_aligns.append(arr[ii:ii + 3])
            else:
                read2_aligns.append(arr[ii:ii + 3])
        # echo(read1_aligns, read2_aligns)

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

    cdef vector[uint64_t] h
    cdef vector[uint16_t] f
    cdef vector[uint32_t] p
    cdef vector[uint16_t] c
    cdef vector[uint64_t] t
    cdef vector[int32_t] cigar_index
    cdef vector[uint32_t] event_pos

    def __cinit__(self):  # Possibly use a vector of structs instead. Reason for doing this way was to save mem
        # node names have the form (mmh3.hash(qname, 42), flag, pos, chrom, tell, cigar index, event pos)
        pass

    cdef void append(self, long a, int b, int c, int d, long e, int f, int g) nogil:
        # echo(a, b, c, d, e, f, g)
        self.h.push_back(a)
        self.f.push_back(b)
        self.p.push_back(c)
        self.c.push_back(d)
        self.t.push_back(e)
        self.cigar_index.push_back(f)
        self.event_pos.push_back(g)

    def __getitem__(self, idx):
        return self.h[idx], self.f[idx], self.p[idx], self.c[idx], self.t[idx], self.cigar_index[idx], self.event_pos[idx]


cdef void update_graph(G, AlignedSegment r, int clip_l, int loci_dist, gettid,
                       overlap_regions, int clustering_dist, clip_scope, PairedEndScoper_t pe_scope,
                       int cigar_index, int event_pos, int paired_end, long tell, genome_scanner,
                       TemplateEdges_t template_edges, NodeToName node_to_name):
                       # unordered_map[int, vector[long_int]]& node_to_name):
    # tt = "HWI-D00360:7:H88WKADXX:1:1106:8344:86803"
    cdef int other_node, clip_left, clip_right
    cdef int current_overlaps_roi, next_overlaps_roi
    cdef bint add_primark_link

    cdef int chrom = r.rname
    cdef int chrom2
    cdef int pos2
    cdef int pos = event_pos #r.pos
    cdef int flag = r.flag
    cdef str qname = r.qname
    cdef uint64_t v

    # if qname == tt:
    #     echo(qname, r.pos, flag, cigar_index, r.cigartuples)
    cdef uint16_t clip_or_wr = 0
    if cigar_index == 0 or cigar_index == len(r.cigartuples) - 1:
        clip_or_wr = 1
    elif 0 < cigar_index < len(r.cigartuples) - 1:
        clip_or_wr = 2

    if paired_end and clip_or_wr <=1 and flag & 8:  # clip event, or whole read, but mate is unmapped
        return
    # Hash qname to save mem
    cdef int node_name = G.addNode()
    ####
    # echo(node_name, flag, pos, chrom, cigar_index, event_pos, qname, clip_or_wr)

    v = xxhasher(bam_get_qname(r._delegate), len(qname), 42)
    # v = mmh3.hash(qname, 42)
    node_to_name.append(v, flag, r.pos, chrom, tell, cigar_index, event_pos)  # Index this list to get the template_name

    genome_scanner.add_to_buffer(r, node_name)  # Add read to buffer


    cdef vector[int] other_nodes

    if paired_end or flag & 1:

        template_edges.add(qname, flag, node_name, r.query_alignment_start)
        # if qname == tt:
        #     echo("template edge", qname, flag, node_name, r.query_alignment_start)
        # Cluster paired-end mates
        pnext = r.pnext
        rnext = r.rnext

        chrom2 = rnext
        pos2 = pnext

        # Special treatment of supplementary and local reads; need to decide where the partner is
        # Either use the rnext:pnext or a site(s) listed in the SA tag: The rnext:pext can often be at the
        # same loci as the read which leads to problems when linking no.2 edges

        add_primark_link = 1

        current_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, r.rname, pos, pos+1)

        if current_overlaps_roi:
            # Cluster soft-clips
            clip_left, clip_right = map_set_utils.clip_sizes(r)
            if clip_left > clip_l:
                best_candidates = clip_scope.update(node_name, r.seq, clip_left, clip_right, chrom, pos)
                if best_candidates is not None:
                    for node_name2 in best_candidates:
                        if not G.hasEdge(node_name, node_name2):
                            G.addEdge(node_name, node_name2, 3)  # weight 3 for black edge
            if clip_right > clip_l:
                best_candidates = clip_scope.update(node_name, r.seq, clip_left, clip_right, chrom, r.reference_end)
                if best_candidates is not None:
                    for node_name2 in best_candidates:
                        if not G.hasEdge(node_name, node_name2):
                            G.addEdge(node_name, node_name2, 3)  # weight 3 for black edge

        if clip_or_wr == 1 and chrom == rnext and abs(pnext - pos) < loci_dist:  # Same loci
            # Use SA tag to generate cluster

            if r.has_tag("SA"):  # Parse SA, first alignment is the other read primary line

                query_length = r.infer_query_length()  # har clips not counted

                for sa_block in r.get_tag("SA").split(";"):  #.split(",", 4)
                    if sa_block == "":
                        break  # End
                    sa = sa_block.split(",", 4)
                    chrom2 = gettid(sa[0])
                    start_pos2 = int(sa[1])
                    # strand = sa[2] == "-"
                    cigar = sa[3]
                    matches = [(int(slen), opp) for slen, opp in re.findall(r'(\d+)([A-Z]{1})', sa[3])]  # parse cigar

                    # Match the current alignment to the SA soft/hard clip
                    diff_left = 10000000000
                    if matches[0][1] in "SH":
                        diff_left = abs(query_length - matches[0][0] - 1)
                    diff_right = 10000000000
                    if matches[len(matches)-1][1] in "SH":
                        diff_right = abs(query_length - matches[len(matches)-1][0] - 1)

                    if diff_left < diff_right:
                        # Use left hand side of SA cigar as break point
                        # if (strand == "-" and flag & 16) or (strand == "+" and not flag & 16):
                        pos2 = start_pos2
                    else:
                        pos2 = start_pos2 + sum(slen for slen, opp in matches if opp == "M" or opp == "D") # todo check inversions

                    add_primark_link = 0
                    next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, chrom2, pos2, pos2+1)

                    if current_overlaps_roi and next_overlaps_roi:
                        continue
                    else:

                        other_nodes = pe_scope.update(node_name, chrom, pos, chrom2, pos2, clip_or_wr)

                        # if tt == r.qname:
                        #     echo("-----", qname, node_name, flag, "a",  pos, pos2, clip_or_wr, other_nodes, cigar_index, r.cigartuples)
                        #     echo(sa)
                        if not other_nodes.empty():
                            for other_node in other_nodes:
                                if not G.hasEdge(node_name, other_node):
                                    G.addEdge(node_name, other_node, 2)
                            other_nodes.clear()


        if add_primark_link: # and paired_end:

            # Use mate information to generate cluster
            if clip_or_wr == 2:  # within read
                chrom2 = chrom
                if r.cigartuples[cigar_index][0] != 1:  # not insertion, use length of cigar event
                    pos2 = event_pos + r.cigartuples[cigar_index][1]
                else:
                    pos2 = event_pos

            elif not clip_or_wr and not flag & 16 and (r.cigartuples[0][0] != 4 or r.cigartuples[0][0] != 5):  # Read on forward strand
                # Use end of read alignment if event likely to be at right side
                pos = r.reference_end

            if current_overlaps_roi and io_funcs.intersecter_int_chrom(overlap_regions, chrom2, pos2, pnext+1):
                # Probably too many reads in ROI to reliably separate out non-soft-clipped reads
                return

            if flag & 2 and not flag & 2048:  # Skip non-discordant/non-split that are same loci
                if abs(r.tlen) < clustering_dist and clip_or_wr != 2:
                    return


            other_nodes = pe_scope.update(node_name, chrom, pos, chrom2, pos2, clip_or_wr)

            # if node_name == 96:
            #     echo("OTHER NODES", other_nodes)
            # if r.cigartuples[cigar_index][0] == 1:
            #     echo("-----", qname, node_name, flag, "b", clip_or_wr,  pos, pos2, other_nodes, cigar_index, r.cigartuples)
            #     quit()
            # if tt == r.qname:
            #     echo("-----", qname, node_name, flag, "b", clip_or_wr,  pos, pos2, other_nodes, cigar_index, r.cigartuples)
                # if r.flag & 64:
                #     quit()
            if not other_nodes.empty():
                for other_node in other_nodes:

                    if not G.hasEdge(node_name, other_node):
                        G.addEdge(node_name, other_node, 2)

    else:  # Single end

        if clip_or_wr == 1:

            template_edges.add(qname, flag, node_name, r.query_alignment_start)  # Add edges between template alignments
            # Use SA tag to get chrom2 and pos2
            if r.has_tag("SA"):  # Parse SA, first alignment for the split read, need to infer the break point for SA

                query_length = r.infer_query_length()  # har clips not counted

                sa = r.get_tag("SA").split(",", 4)
                chrom2 = gettid(sa[0])
                start_pos2 = int(sa[1])
                strand = sa[2] == "-"
                cigar = sa[3]
                matches = [(int(slen), opp) for slen, opp in re.findall(r'(\d+)([A-Z]{1})', sa[3])]  # parse cigar

                # Match the current alignment to the SA soft/hard clip
                diff_left = 10000000000
                if matches[0][1] in "SH":
                    diff_left = abs(query_length - matches[0][0] - 1)
                diff_right = 10000000000
                if matches[len(matches)-1][1] in "SH":
                    diff_right = abs(query_length - matches[len(matches)-1][0] - 1)

                if diff_left < diff_right:
                    # Use left hand side of SA cigar as break point
                    # if (strand == "-" and flag & 16) or (strand == "+" and not flag & 16):
                    pos2 = start_pos2
                else:
                    pos2 = start_pos2 + sum(slen for slen, opp in matches if opp == "M" or opp == "D") # todo check for inversions

                # Find position on reference
                other_nodes = pe_scope.update(node_name, chrom, event_pos, chrom2, pos2, clip_or_wr)
                if not other_nodes.empty():
                    for other_node in other_nodes:
                        if not G.hasEdge(node_name, other_node):
                            # echo("sa edge", node_name, other_node, cigar_index, len(r.cigartuples), qname)
                            G.addEdge(node_name, other_node, 2)
        elif clip_or_wr == 2:  # Sv within read
            chrom2 = r.rname
            if r.cigartuples[cigar_index][0] != 1:
                pos2 = event_pos + r.cigartuples[cigar_index][1]
            else:
                pos2 = event_pos

            # if node_name == 1:
            #     echo(chrom, event_pos, chrom2, pos2,)
            other_nodes = pe_scope.update(node_name, chrom, event_pos, chrom2, pos2, clip_or_wr)
            if not other_nodes.empty():
                for other_node in other_nodes:
                    if not G.hasEdge(node_name, other_node):
                        # echo("wr edge", node_name, other_node)
                        G.addEdge(node_name, other_node, 2)


cpdef tuple construct_graph(genome_scanner, infile, int max_dist, int clustering_dist, int k=16, int m=7, int clip_l=21,
                            int minimizer_support_thresh=2, int minimizer_breadth=3,
                            int minimizer_dist=10, int mapq_thresh=1, debug=None, int min_support=3, procs=1,
                            int paired_end=1):

    t0 = time.time()
    click.echo("Building graph with clustering distance {} bp, scope length {} bp".format(max_dist, clustering_dist),
               err=True)

    cdef TemplateEdges_t template_edges = TemplateEdges()  # Edges are added between alignments from same template, after building main graph
    # cdef unordered_map[int, vector[long_int]] node_to_name
    cdef int event_pos, cigar_index, opp, length, added

    node_to_name = NodeToName()  # Map of nodes -> read ids

    clip_scope = ClipScoper(minimizer_dist, k=k, m=m, clip_length=clip_l,  # Keeps track of local reads
                       minimizer_support_thresh=minimizer_support_thresh,
                       minimizer_breadth=minimizer_breadth)

    # Infers long-range connections, outside local scope using pe information
    cdef PairedEndScoper_t pe_scope = PairedEndScoper(max_dist, clustering_dist, infile.header.nreferences)
    cdef int loci_dist = max_dist #int(max_dist * 1.5)
    cdef long tell
    # cdef Py_SimpleGraph_t G = map_set_utils.Py_SimpleGraph()
    G = map_set_utils.Py_SimpleGraph()
    overlap_regions = genome_scanner.overlap_regions  # Get overlapper, intersect reads with intervals

    gettid = infile.gettid

    # debug = "HISEQ1:12:H8GVUADXX:2:2109:12037:62696"

    debug_nodes = set([])

    for chunk in genome_scanner.iter_genome():

        for r, tell in chunk:

            if r.mapq < mapq_thresh:
                continue
            # echo(r.qname, r.cigartuples, r.has_tag("SA"))
            event_pos = r.pos

            added = 0
            if len(r.cigartuples) > 1:

                for cigar_index, (opp, length) in enumerate(r.cigartuples):

                    if (opp == 4 or opp == 5) and (length > 30 or r.has_tag("SA")):  #length > 30:
                        update_graph(G, r, clip_l, loci_dist, gettid,
                           overlap_regions, clustering_dist, clip_scope, pe_scope,
                           cigar_index, event_pos, paired_end, tell, genome_scanner,
                           template_edges, node_to_name)
                        added += 1

                    elif opp == 1 and length > 30:
                        update_graph(G, r, clip_l, loci_dist, gettid,
                           overlap_regions, clustering_dist, clip_scope, pe_scope,
                           cigar_index, event_pos, paired_end, tell, genome_scanner,
                           template_edges, node_to_name)
                        added += 1

                    elif opp == 2:
                        if length > 30:
                            update_graph(G, r, clip_l, loci_dist, gettid,
                               overlap_regions, clustering_dist, clip_scope, pe_scope,
                               cigar_index, event_pos, paired_end, tell, genome_scanner,
                               template_edges, node_to_name)

                            added += 1
                        event_pos += length

                    else:
                        event_pos += length

            if not added:
                # -1 means whole alignment
                cigar_index = -1
                update_graph(G, r, clip_l, loci_dist, gettid,
                           overlap_regions, clustering_dist, clip_scope, pe_scope,
                           cigar_index, event_pos, paired_end, tell, genome_scanner,
                           template_edges, node_to_name)


    add_template_edges(G, template_edges)

    return G, node_to_name


cpdef tuple get_block_components(G, node_to_name, infile, read_buffer,
                         int min_support):
                         #debug_nodes):

    # Turn graph into a block model, first split into connected components,
    # Then for each component, split into block nodes which are locally interacting nodes (black and grey edges)
    # block nodes edges result from collapsing white edges

    cdef dict partitions, support_within, reads
    cdef int v, item_idx, item

    t0 = time.time()
    cmp = G.connectedComponents()  # Flat vector, components are separated by -1

    cc = []
    last_i = 0
    for item_idx, item in enumerate(cmp):
        if item == -1:
            cc.append((last_i, item_idx))
            last_i = item_idx + 1

    return cmp, cc


cpdef dict get_reads(infile, sub_graph_reads):

    rd = dict()
    # coords = []
    # missing = set([])  # Integer nodes
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
        # n1 = (mmh3.hash(a.qname, 42), a.flag, a.pos, a.rname, p)
        n1 = (v, a.flag, a.pos, a.rname, p)
        # Try next few reads, sometimes they are on top of one another
        if n1 != node:
            for j in range(5):
                a = next(infile)
                n2 = (xxhasher(bam_get_qname(a._delegate), len(a.qname), 42), a.flag, a.pos, a.rname, p)
                if n2 == node:
                    rd[int_node] = a
                    break
            # else:
            #     missing.add(int_node)

        else:
            rd[int_node] = a

    return rd


cdef set BFS_local(G, int source, unordered_set[int]& visited ):

    # Mark all the vertices as not visited
    # visited = set([])

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

        # visited.add(u)
        visited.insert(u)

    return nodes_found #, visited


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
        # for v in parts[0]:
        #     echo(f"--->node={v}", len(list(G.neighbors(int(v)))))
        #     for v2 in G.neighbors(int(v)):
        #
        #         echo(f"node={v}", f"node2={v2}", f"w={G.weight(int(v), int(v2))}")

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

                # j = p2i[child]  # Partition of neighbor node
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

    cdef first, second
    for t in seen_t:
        if sum(len(item) for item in counts[t]) < min_support:  # .values()
            del counts[t]
            first = t[0]
            second = t[1]
            if first in self_counts and len(self_counts[first]) < min_support:
                del self_counts[first]
            if second in self_counts and len(self_counts[second]) < min_support:
                del self_counts[second]

    return counts, self_counts


cpdef dict proc_component(node_to_name, component, read_buffer, infile, G,
                         int min_support):

    n2n = {}
    reads = {}

    cdef int v
    for v in component:
        if v in read_buffer:
            reads[v] = read_buffer[v]
        # Need to keep a record of all node info, and cigar indexes
        n2n[v] = node_to_name[v]

    # Explore component for locally interacting nodes; create partitions using these
    partitions = get_partitions(G, component)
    support_between, support_within = count_support_between(G, partitions, min_support)

    if len(support_between) == 0 and len(support_within) == 0:
        if len(reads) >= min_support:
            # echo("1parts", partitions)
            # echo("1s_between", support_between)
            # echo("1s_within", support_within)
            # echo("1n2n", n2n.keys())
            return {"parts": {}, "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}
        else:
            return {}

    sb = {}
    kept = set([])
    deleted = set([])
    for edge, vd in support_between.items():

        sup = sum([len(vv) for vv in vd])  # .values()
        if sup >= min_support:
            sb[edge] = vd
            kept.add(edge[0])
            kept.add(edge[1])
        else:
            deleted.add(edge[0])
            deleted.add(edge[1])

    deleted -= kept

    # check deleted for support within, discard if min_support not reached
    for block_node in deleted:
        if block_node in support_within:
            if len(support_within[block_node]) < min_support:
                del partitions[block_node]
        else:
            del partitions[block_node]

    # echo("parts", partitions)
    # echo("s_between", sb)
    # echo("s_within", support_within)
    # echo("n2n", n2n.keys())

    return {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}

