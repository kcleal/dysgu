#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

from __future__ import absolute_import
import time
import click
from collections import defaultdict, deque
import mmh3

import sortedcontainers
import ncls

import pysam
import numpy as np
cimport numpy as np
from cpython cimport array
import array

from redblackpy.series.tree_series import Series
from redblackpy.series.tree_series cimport Series

from libcpp.vector cimport vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string

from libc.stdint cimport uint32_t

from cython.operator import dereference, postincrement

# from libcpp.set cimport set as cpp_set
# from libcpp.unordered_set cimport unordered_set as cpp_u_set

# from libcpp.string cimport string as cpp_string
# from libc.stdint cimport int32_t, int64_t
# from cython.operator cimport dereference as deref, preincrement as inc

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK
# from libc.stdint cimport uint32_t, uint8_t

from dysgu cimport map_set_utils
from dysgu import io_funcs

from dysgu.map_set_utils cimport robin_set


ctypedef cpp_pair[int, int] cpp_item
ctypedef map_set_utils.Py_IntSet Py_IntSet
ctypedef map_set_utils.Py_Int2IntMap Py_Int2IntMap

ctypedef map_set_utils.Py_SimpleGraph Py_SimpleGraph

ctypedef cpp_map[int, cpp_item] ScopeItem_t
ctypedef vector[ScopeItem_t] ForwardScope_t


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
        cdef int i
        with nogil:
            for i in range(a.shape[0]):
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

    # cdef cpp_set[int] seen2  # Using
    # cdef cpp_u_set[int] seen2
    seen2 = set([])

    # cdef cpp_item last
    for i in range(end):

        hx2 = mmh3.hash(s[i:i+m], 42)

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


cdef struct LocalItem:
    int pos
    int chrom2
    int pos2
    int node_name


cdef LocalItem make_local_item(int pos, int chrom2, int pos2, int node_name):
    cdef LocalItem item
    item.pos = pos
    item.chrom2 = chrom2
    item.pos2 = pos2
    item.node_name = node_name
    return item


# cdef class SortItem:
#     cdef readonly int pos2, node_name, current_pos
#     def __cinit__(self, pos2, node_name, current_pos):
#         self.pos2 = pos2
#         self.node_name = node_name
#         self.current_pos = current_pos
#     def __lt__(self, other):
#         return self.pos2 < other
#     def __gt__(self, other):
#         return other <= self.pos2
#     def __eq__(self, other):
#         return self.pos2 == other
#     def __repr__(self):
#         return f"({self.pos2}, {self.node_name}, {self.current_pos})"
#
#
# cdef int check_distance(SortItem v, int current_pos, int pos2, int max_dist):
#     if abs(v.current_pos - current_pos) < max_dist and abs(v.pos2 - pos2) < max_dist:
#         return 1
#     return 0



cdef class PairedEndScoper:

    cdef float clst_dist
    cdef int max_dist
    cdef int local_chrom
    # cdef list forward_scope
    cdef list forward_scope2

    cdef cpp_deque[LocalItem] local_scope  # in current region (pos, chrom2)

    def __init__(self, max_dist, n_references):
        # Keep reads within this distance in scope, *5 seems to work well but not sure why
        self.clst_dist = max_dist * 5
        self.max_dist = max_dist
        self.local_chrom = -1

        # self.forward_scope = [sortedcontainers.SortedList() for i in range(n_references)]
        self.forward_scope2 = [Series(dtype="object") for i in range(n_references)]

    cdef int check_far(self, int v_current, int v_pos2, int current_pos, int pos2):
        if abs(v_current - current_pos) < self.max_dist and abs(v_pos2 - pos2) < self.max_dist:
            return 1
        return 0

    cdef int check_close(self, int v_current, int current_pos):
        if abs(v_current - current_pos) < self.max_dist:
            return 1
        return 0

    cdef list update(self, int node_name, int current_chrom, int current_pos, int chrom2, int pos2):

        if current_chrom != self.local_chrom:

            self.local_chrom = current_chrom

            self.local_scope.clear()
            self.local_scope.push_back(make_local_item(current_pos, chrom2, pos2, node_name))

            # for idx in range(len(self.forward_scope)):
            #     if len(self.forward_scope[idx]) > 0:
            #         self.forward_scope[idx] = sortedcontainers.SortedList()
            # self.forward_scope[chrom2].add(SortItem(pos2, node_name, current_pos))

            for idx in range(len(self.forward_scope2)):
                if len(self.forward_scope2[idx]) > 0:
                    self.forward_scope2[idx] = Series(dtype="object")
            self.forward_scope2[chrom2][pos2] = (node_name, current_pos)
            return []

        cdef LocalItem item
        while self.local_scope.size() > 0:

            item = self.local_scope[0]
            if abs(item.pos - current_pos) > self.clst_dist * 5:

                # self.forward_scope[item.chrom2].remove(item.pos2)

                self.forward_scope2[item.chrom2].erase(item.pos2)

                self.local_scope.pop_front()
            else:
                break

        # TODO, use multiple values at each pos; use list instead
        #
        found2 = []
        t = self.forward_scope2[chrom2]
        len_t = len(t)
        if len_t == 0:
            t.insert(pos2, (node_name, current_pos))
            return []


        elif len_t == 1:
            v = t.ceil(pos2)

            if v[0] == pos2:
                self.forward_scope2[chrom2][pos2] = (node_name, current_pos)
                if self.check_close(v[1][1], current_pos):
                    found2.append(v[1][0])

            else:
                if v[0] is None:
                    v = t.floor(pos2 - 1)

                if self.check_far(v[1][1], v[0], current_pos, pos2):
                    found2.append(v[1][0])
                self.forward_scope2[chrom2].insert(pos2, (node_name, current_pos))

        else:
            lower = t.floor(pos2 - 1)
            upper = t.ceil(pos2)

            if upper[0] == pos2:

                if self.check_far(upper[1][1], upper[0], current_pos, pos2):
                    found2.append(upper[1][0])

                # Key exists, test upper and lower again
                upper2 = t.ceil(pos2 + 1)

                if upper2[0] is not None and self.check_close(upper2[1][1], current_pos):
                    found2.append(upper2[1][0])
                if lower[0] is not None and self.check_close(lower[1][1], current_pos):
                    found2.append(lower[1][0])
                # update existing key
                self.forward_scope2[chrom2][pos2] = (node_name, current_pos)

            else:
                if upper[0] is None:  # look for next closest
                    upper = t.floor(lower[0] - 1)

                elif lower[0] is None:
                    lower = t.ceil(upper[0] + 1)

                if self.check_far(lower[1][1], lower[0], current_pos, pos2):
                    found2.append(lower[1][0])
                if self.check_far(upper[1][1], upper[0], current_pos, pos2):
                    found2.append(upper[1][0])

                self.forward_scope2[chrom2].insert(pos2, (node_name, current_pos))

        # Find other nodes in forward_scope
        # t = self.forward_scope[chrom2]
        # len_t = len(t)
        # if len_t == 0:
        #     t.add(SortItem(pos2, node_name, current_pos))
        #     return []
        #
        # idx = t.bisect_left(pos2)
        # found = []
        # if idx == 0:
        #     if check_distance(t[idx], current_pos, pos2, self.max_dist):
        #         found.append(t[idx].node_name)
        #     if len(t) > 1:
        #         if check_distance(t[1], current_pos, pos2, self.max_dist):
        #             found.append(t[1].node_name)
        #
        #
        # elif idx == len_t:
        #     if check_distance(t[idx - 1], current_pos, pos2, self.max_dist):
        #         found.append(t[idx - 1].node_name)
        #     if len(t) > 1:
        #         if check_distance(t[idx - 2], current_pos, pos2, self.max_dist):
        #             found.append(t[idx - 2].node_name)
        #
        # else:  # Add upstream and downstream
        #     if check_distance(t[idx], current_pos, pos2, self.max_dist):
        #         found.append(t[idx].node_name)
        #
        #     if check_distance(t[idx - 1], current_pos, pos2, self.max_dist):
        #         found.append(t[idx - 1].node_name)
        #
        # self.forward_scope[chrom2].add(SortItem(pos2, node_name, current_pos))


        self.local_scope.push_back(make_local_item(current_pos, chrom2, pos2, node_name))


        # echo(found, found2)
        # echo("pos2", pos2, "current_pos", current_pos)
        # echo("ori", self.forward_scope[chrom2])
        # echo("new", "len", len(self.forward_scope2[chrom2]), self.forward_scope2[chrom2])
        #
        # if len(found2) == 2 and found2[0] == found2[1]:
        #     quit()
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



cdef void add_template_edges(G, TemplateEdges template_edges): #TemplateEdges template_edges):

    cdef int ii, u_start, v_start, u, v, uflag, vflag

    # for qname, arr in template_edges.templates.items():  # normally 2 reads, or >2 if supplementary reads
    # for qname, (read1_aligns, read2_aligns) in template_edges.templates.items():
    for qname, arr in template_edges.iterate_map():
    #     echo(arr)
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

    cdef array.array h
    cdef array.array f
    cdef array.array p
    cdef array.array c
    cdef array.array t

    def __cinit__(self):  # Possibly use a vector of structs instead. Reason for doing this way was to save mem
        # node names have the form (mmh3.hash(qname, 42), flag, pos, chrom, tell)
        self.h = array.array("l", [])  # signed long
        self.f = array.array("H", [])  # unsigned short
        self.p = array.array("L", [])  # unsigned long
        self.c = array.array("H", [])
        self.t = array.array("l", [])

    cdef void append(self, a, b, c, d, e):
        self.h.append(a)
        self.f.append(b)
        self.p.append(c)
        self.c.append(d)
        self.t.append(e)

    def __getitem__(self, idx):
        return self.h[idx], self.f[idx], self.p[idx], self.c[idx], self.t[idx]


def construct_graph(genome_scanner, infile, int max_dist, int clustering_dist, int k=16, int m=7, int clip_l=21,
                    int minimizer_support_thresh=2, int minimizer_breadth=3,
                    int minimizer_dist=10, int mapq_thresh=10, debug=None, int min_support=3, procs=1,
                    int paired_end=1):

    t0 = time.time()
    click.echo("Building graph with clustering distance {} bp".format(max_dist), err=True)

    template_edges = TemplateEdges()  # Edges are added between alignments from same template, after building main graph

    node_to_name = NodeToName()  # Map of useful nodes -> read names

    clip_scope = ClipScoper(minimizer_dist, k=k, m=m, clip_length=clip_l,  # Keeps track of local reads
                       minimizer_support_thresh=minimizer_support_thresh,
                       minimizer_breadth=minimizer_breadth)


    pe_scope = PairedEndScoper(max_dist, infile.header.nreferences)  # Infers long-range connections, outside local scope using pe information

    overlap_regions = genome_scanner.overlap_regions  # Get overlapper, intersect reads with intervals

    gettid = infile.gettid

    cdef int flag, pos, rname, pnext, rnext, node_name, node_name2, chrom, chrom2, clip_left, clip_right
    cdef str qname

    cdef int ol_include, add_primary_link
    cdef int current_overlaps_roi, next_overlaps_roi

    cdef int loci_dist = int(max_dist * 1.5)
    t0 = time.time()

    cdef Py_SimpleGraph G = map_set_utils.Py_SimpleGraph()

    # debug = "HWI-D00360:5:H814YADXX:2:1110:10283:75050"

    debug_nodes = set([])
    for chunk in genome_scanner.iter_genome():

        for r, tell in chunk:

            if r.mapq < mapq_thresh:
                continue

            chrom = r.rname
            pos = r.pos
            flag = r.flag
            qname = r.qname

            # Hash qname to save mem
            node_name = G.addNode()
            node_to_name.append(mmh3.hash(qname, 42), flag, pos, chrom, tell)  # Index this list to get the template_name

            genome_scanner.add_to_buffer(r, node_name)  # Add read to buffer
            template_edges.add(qname, flag, node_name, r.query_alignment_start)

            if paired_end:

                # Cluster paired-end mates
                pnext = r.pnext
                rnext = r.rnext

                # Special treatment of supplementary and local reads; need to decide where the partner is
                # Either use the rnext:pnext or a site(s) listed in the SA tag: The rnext:pext can often be at the
                # same loci as the read which leads to problems when linking no.2 edges

                add_primark_link = 1

                current_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, r.rname, pos, pos+1)

                if current_overlaps_roi:
                    # Cluster soft-clips
                    clip_left, clip_right = map_set_utils.clip_sizes(r)
                    if clip_left > clip_l or clip_right > clip_l:
                        best_candidates = clip_scope.update(node_name, r.seq, clip_left, clip_right, chrom, pos)

                        if best_candidates is not None:
                            for node_name2 in best_candidates:
                                if not G.hasEdge(node_name, node_name2):
                                    G.addEdge(node_name, node_name2, 3)  # weight 3 for black edge

                if chrom == rnext and abs(pnext - pos) < loci_dist:  # Same loci
                    # Use SA tag to generate cluster
                    if r.has_tag("SA"):  # Parse SA, first alignment is the other read primary line
                        sa = r.get_tag("SA").split(",", 2)
                        chrom2 = gettid(sa[0])
                        pos2 = int(sa[1])

                        add_primark_link = 0
                        next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, chrom2, pos2, pos2+1)

                        if current_overlaps_roi and next_overlaps_roi:
                            continue
                        else:
                            for other_node in pe_scope.update(node_name, chrom, pos, chrom2, pos2):
                                if not G.hasEdge(node_name, other_node):
                                    G.addEdge(node_name, other_node, 2)


                if add_primark_link == 1: # and paired_end:
                    # Use mate information to generate cluster
                    next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, rnext, pnext, pnext+1)

                    if current_overlaps_roi and next_overlaps_roi:
                        # Probably too many reads in ROI to reliably separate out non-soft-clipped reads
                        continue

                    if flag & 2 and not flag & 2048:  # Skip non-discordant that are same loci
                        if abs(r.tlen) < clustering_dist:
                           continue

                    other_nodes = pe_scope.update(node_name, chrom, pos, rnext, pnext)

                    if other_nodes:
                        for other_node in other_nodes:
                            # if node_to_name[other_node] not in templates_connected and \
                            if not G.hasEdge(node_name, other_node):
                                G.addEdge(node_name, other_node, 2)

            else:  # Single end
                # Use SA tag to generate cluster
                if r.has_tag("SA"):  # Parse SA, first alignment is the primary?
                    sa = r.get_tag("SA").split(",", 2)
                    chrom2 = gettid(sa[0])
                    pos2 = int(sa[1])

                    for other_node in pe_scope.update(node_name, chrom, pos, chrom2, pos2):
                        if not G.hasEdge(node_name, other_node):

                            G.addEdge(node_name, other_node, 2)

                else:
                    pass
                    # Check for SVs within reads from cigar
                    # current_pos = r.pos
                    # for opp, length in r.cigartuples:
                    #     if (opp == 2 or opp == 1) and length > 50:
                    #         if opp == 2:
                    #             pos2 = current_pos + length + 1
                    #
                    #             for other_node in pe_scope.update(node_name, chrom, pos, chrom2, pos2):
                    #                 if not G.hasEdge(node_name, other_node):
                    #                     G.addEdge(node_name, other_node, 2)
                    #
                    #             current_pos = pos2
                    #
                    #     if opp == 1:
                    #         current_pos += 1
                    #
                    #     elif opp == 0 or opp == 7 or opp == 8 or opp == 3:
                    #         current_pos += length



    click.echo(f"Processed minimizers {time.time() - t0}", err=True)

    t2 = time.time()

    add_template_edges(G, template_edges)

    # Free mem
    # template_edges = None
    # pe_scope = None
    click.echo(f"Paired end info in {time.time() - t2}", err=True)

    return G, node_to_name



    # cdef dict component



    # read_buffer = genome_scanner.read_buffer
    #
    # click.echo(f"Constructed graph, edges={G.edgeCount()}, {int(time.time() - t0)} s, processing block model", err=True)
    # for component in get_block_components(G, node_to_name, infile, read_buffer, min_support):
    #                                       # debug_nodes):
    #     yield component



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


    # cdef int cnt = 0
    # for start_i, end_i in cc:
    #     cnt += 1
    #     if end_i - start_i >= min_support:
    #         component = list(cmp[start_i: end_i])
    #
    #
    #         res = proc_component(node_to_name, component, read_buffer, infile, G, min_support)
    #
    #         if res:
    #             yield res
    #     # Reduce size of graph
    #     # for v in component:
    #     #     G.removeNode(v)
    #
    # click.echo(f"Processed components n={cnt} {int(time.time() - t0)} s", err=True)


cpdef dict get_reads(infile, sub_graph_reads):

    rd = dict()
    # coords = []
    # missing = set([])  # Integer nodes
    cdef int j, int_node
    cdef long int p

    for int_node, node in sub_graph_reads.items():

        p = node[4]
        infile.seek(p)
        a = next(infile)
        n1 = (mmh3.hash(a.qname, 42), a.flag, a.pos, a.rname, p)

        # Try next few reads, sometimes they are on top of one another
        if n1 != node:
            for j in range(5):
                a = next(infile)
                n2 = (mmh3.hash(a.qname, 42), a.flag, a.pos, a.rname, p)
                if n2 == node:
                    rd[int_node] = a
                    break
            # else:
            #     missing.add(int_node)

        else:
            rd[int_node] = a

    return rd


cdef set BFS_local(G, int source, robin_set[int]& visited ):

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

    cdef robin_set[int] seen

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
    # p2i = dict()
    cdef Py_Int2IntMap p2i = map_set_utils.Py_Int2IntMap()
    for i, p in parts.items():
        for node in p:
            p2i.insert(node, i)

        # p2i.update({node: i for node in p})

    # Count the links between partitions. Split reads into sets ready for calling
    # No counting of read-pairs templates or 'support', just a count of linking alignments
    # counts (part_a, part_b): {part_a: {node 1, node 2 ..}, part_b: {node4, ..} }
    counts = {}

    self_counts = {}  # defaultdict(empty_long_array)

    # cdef Py_IntSet seen_nodes = map_set_utils.Py_IntSet()

    seen_t = set([])
    for i, p in parts.items():

        current_t = set([])
        for node in p:
            any_out_edges = 0  # Keeps track of number of outgoing pairs, or self edges

            # if seen_nodes.has_key(node):
            #     continue

            for child in G.neighbors(node):
                if not p2i.has_key(child): # or seen_nodes.has_key(child):  # child in p2i
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

                    # counts[t][i].add(node)
                    # counts[t][j].add(child)

                    current_t.add(t)

            # seen_nodes.insert(node)

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
        else:
            n2n[v] = node_to_name[v]

    # Explore component for locally interacting nodes; create partitions using these
    partitions = get_partitions(G, component)
    support_between, support_within = count_support_between(G, partitions, min_support)

    if len(support_between) == 0 and len(support_within) == 0 and len(reads) >= min_support:
        return {"parts": {}, "s_between": {}, "reads": reads, "s_within": {}, "n2n": n2n}

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

    return {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}

