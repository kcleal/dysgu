#cython: language_level=3, boundscheck=False, c_string_type=unicode, c_string_encoding=utf8, infer_types=True

from __future__ import absolute_import
import time
import click
from collections import defaultdict, deque

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
from libc.math cimport fabs as c_fabs

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


cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil:
    # Insertions have same x1/y1 position
    if x1 == x2 or y1 == y2:
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
    if (overlap / float(c_abs(x2 - x1))) > 0.1 and (overlap / float(c_abs(y2 - y1))) > 0.1:
        return True


cdef bint span_position_distance(int x1, int x2, int y1, int y2) nogil:
    # https://github.com/eldariont/svim/blob/master/src/svim/SVIM_clustering.py
    cdef int span1, span2, max_span
    cdef float span_distance, position_distance, center1, center2
    if x1 == x2:
        span1 = 1
        center1 = x1
    else:
        span1 = c_abs(x2 - x1)
        center1 = (x1 + x2) / 2
    if y1 == y2:
        span2 = 1
        center2 = y2
    else:
        span2 = c_abs(y2 - y1)
        center2 = (y1 + y2) / 2
    max_span = max(span1, span2)
    position_distance = c_fabs(center1 - center2) # 1 #distance_normalizer
    span_distance = <float>c_abs(span1 - span2) / max_span
    # echo("pd", position_distance, center1, center2)
    # echo((position_distance / max_span), span_distance, center1, center2 )
    if (position_distance / max_span) < 0.2 and span_distance < 0.3:
    # if position_distance < 100 and span_distance < 0.08:
        return 1
    return 0


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

    cdef vector[int] find_other_nodes(self, int node_name, int current_chrom, int current_pos, int chrom2, int pos2,
                            int clip_or_wr) nogil:

        cdef int idx, i, count_back, steps, node_name2 #, lowest_pos
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

            # local_it = forward_scope.begin()
            # while local_it != forward_scope.end():
            #     vitem = dereference(local_it)
            #     echo("before", vitem.first, vitem.second)
            #     preincrement(local_it)


            # lowest_pos = current_pos
            # if current_chrom == chrom2 and pos2 < current_pos:  # Make sure only clear up to lowest position
            #     lowest_pos = pos2

            local_it = self.loci.lower_bound(current_pos - self.clst_dist)
            if local_it != self.loci.begin():
                # echo("dropping up to", dereference(local_it).first)
                self.loci.erase(self.loci.begin(), local_it)

            if current_chrom != chrom2:  # Erase items out of range in forward scope
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

                    if (clip_or_wr == 2 and vitem.second.clip_or_wr == 3) or (clip_or_wr == 3 and vitem.second.clip_or_wr == 2):
                        steps += 1
                        continue  # Dont connect insertions to deletions

                    node_name2 = vitem.second.node_name
                    if node_name2 != node_name:  # Can happen due to within-read events
                        # echo("1", current_pos, pos2, is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2), span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2))
                        if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):

                            sep = c_abs(vitem.first - pos2)
                            sep2 = c_abs(vitem.second.pos2 - current_pos)
                            # echo(sep, sep2, self.max_dist, vitem.second.clip_or_wr)
                            if sep < self.max_dist and sep2 < self.max_dist:
                                if sep < 35 and (clip_or_wr > 0 or vitem.second.clip_or_wr):
                                    found_exact.push_back(node_name2)
                                else:
                                    found2.push_back(node_name2)
                            elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2):
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
                        if (clip_or_wr == 2 and vitem.second.clip_or_wr == 3) or (clip_or_wr == 3 and vitem.second.clip_or_wr == 2):
                            steps += 1
                            continue  # Dont connect insertions to deletions
                        node_name2 = vitem.second.node_name
                        if node_name2 != node_name:
                            # echo("2", is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2))
                            if current_chrom != chrom2 or is_reciprocal_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2): # is_overlapping(current_pos, pos2, vitem.first, vitem.second.pos2):
                                sep = c_abs(vitem.first - pos2)
                                sep2 = c_abs(vitem.second.pos2 - current_pos)

                                if sep < self.max_dist and vitem.second.chrom2 == chrom2 and \
                                        sep2 < self.max_dist:
                                    if sep < 35 and (clip_or_wr > 0 or vitem.second.clip_or_wr):
                                        found_exact.push_back(node_name2)
                                    else:
                                        found2.push_back(node_name2)
                                elif span_position_distance(current_pos, pos2, vitem.first, vitem.second.pos2):
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
                            int clip_or_wr) nogil:

        # Add to scopes, if event is within read, add two references to forward scope. Otherwise when the
        # first break point drops from the local scope, the forward scope break will not connect with other
        # events that are added later. This is a fix for long read deletions mainly
        cdef cpp_map[int, LocalVal]* forward_scope = &self.chrom_scope[chrom2]
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

        if clip_or_wr == 2:  # Deletion type
            local_it = forward_scope.find(current_pos)
            if local_it == forward_scope.end():
                dereference(forward_scope)[current_pos] = make_local_val(chrom2, pos2, node_name, clip_or_wr)
            else:
                dereference(local_it).second = make_local_val(chrom2, pos2, node_name, clip_or_wr)


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
    # Index these vectors to get the unique 'template_name'
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
        self.h.push_back(a)
        self.f.push_back(b)
        self.p.push_back(c)
        self.c.push_back(d)
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
    return start, end, r.pos, ref_end, r.rname, 0


cdef alignments_from_sa_tag(r, gettid):
    # Puts other alignments in order of query position, gets index of the current alignment
    cdef int query_length = r.infer_read_length()  # Note, this also counts hard-clips
    cdef int chrom2, start_pos2, query_start, query_end, ref_start, ref_end, start_temp

    current_strand = "-" if r.flag & 16 else "+"
    query_aligns = [get_query_pos_from_cigartuples(r)]

    for sa_block in r.get_tag("SA").split(";"):
        if sa_block == "":
            break  # End
        sa = sa_block.split(",", 4)
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

        query_aligns.append((query_start, query_end, ref_start, ref_end, chrom2))

    query_aligns = sorted(query_aligns)
    cdef int index = 0
    for index in range(len(query_aligns)):
        if len(query_aligns[index]) == 6:
            break
    # if r.qname == "m64004_190803_004451/6359627/ccs":
    #     echo(query_aligns, index)
    #     quit()
    return query_aligns, index


cdef connect_right(a, b):
    event_pos = a[3]
    chrom = a[4]
    pos2 = b[2]
    chrom2 = b[4]
    return event_pos, chrom, pos2, chrom2


cdef connect_left(a, b):
    event_pos = a[2]
    chrom = a[4]
    pos2 = b[3]
    chrom2 = b[4]
    return event_pos, chrom, pos2, chrom2


cdef bint add_to_graph(G, AlignedSegment r, PairedEndScoper_t pe_scope, TemplateEdges_t template_edges,
                       NodeToName node_to_name, genome_scanner,
                       int flag, int chrom, long tell, int cigar_index, int event_pos,
                       int chrom2, int pos2, int clip_or_wr):

    # Adds relevant information to graph and other data structures for further processing

    cdef vector[int] other_nodes  # Other alignments to add edges between
    cdef int node_name = G.addNode()
    cdef uint64_t v = xxhasher(bam_get_qname(r._delegate), len(r.qname), 42)  # Hash qname to save mem

    node_to_name.append(v, flag, r.pos, chrom, tell, cigar_index, event_pos)  # Index this list to get the template_name

    genome_scanner.add_to_buffer(r, node_name)  # Add read to buffer

    if clip_or_wr < 2:  # Prevents joining up within-read svs within between-read svs
        template_edges.add(r.qname, flag, node_name, r.query_alignment_start)

    other_nodes = pe_scope.find_other_nodes(node_name, chrom, event_pos, chrom2, pos2, clip_or_wr)

    pe_scope.add_item(node_name, chrom, event_pos, chrom2, pos2, clip_or_wr)

    # echo(r.qname, node_name, event_pos, pos2)
    # if r.qname == "m64004_190803_004451/6359627/ccs":
    #       echo("@", node_name, event_pos, pos2, list(other_nodes), cigar_index)

    if not other_nodes.empty():
        for other_node in other_nodes:
            if not G.hasEdge(node_name, other_node):
                # 2 signifies that this is a local edge as opposed to a template edge
                G.addEdge(node_name, other_node, 2)
        return 1

    return 0


cdef void process_alignment(G, AlignedSegment r, int clip_l, int loci_dist, gettid,
                       overlap_regions, int clustering_dist, PairedEndScoper_t pe_scope,
                       int cigar_index, int event_pos, int paired_end, long tell, genome_scanner,
                       TemplateEdges_t template_edges, NodeToName node_to_name,
                            int clip_or_wr, int cigar_pos2):

    # Determines where the break point on the alignment is before adding to the graph
    cdef int other_node, clip_left, clip_right
    cdef int current_overlaps_roi, next_overlaps_roi
    cdef bint add_primark_link
    cdef int chrom = r.rname
    cdef int chrom2
    cdef int pos2
    cdef int flag = r.flag
    cdef str qname = r.qname
    cdef uint64_t v

    if paired_end and clip_or_wr == 1 and flag & 8:  # clip event, or whole read, but mate is unmapped
        return

    cdef bint success

    if paired_end or flag & 1:

        # Cluster paired-end mates
        pnext = r.pnext
        rnext = r.rnext
        chrom2 = rnext
        pos2 = pnext

        # Special treatment of supplementary and local reads; need to decide where the partner is
        # Either use the rnext:pnext or a site(s) listed in the SA tag: The rnext:pext can often be at the
        # same loci as the read which leads to problems when linking w=2 edges

        add_primary_link = 1
        current_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, r.rname, event_pos, event_pos+1)

        if clip_or_wr == 1 and chrom == rnext and abs(pnext - event_pos) < loci_dist:  # Same loci

            # Parse SA tag. For paired reads
            if r.has_tag("SA"):  # Parse SA, first alignment is the other read primary line

                all_aligns, index = alignments_from_sa_tag(r, gettid)
                event = all_aligns[index]

                if len(all_aligns) == 1:
                    return  # shouldnt happen

                add_primary_link = 0

                if index < len(all_aligns) - 1:  # connect to next
                    event_pos, chrom, pos2, chrom2 = connect_right(all_aligns[index], all_aligns[index + 1])
                    cigar_index = len(r.cigartuples) - 1
                    success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                           tell, cigar_index, event_pos, chrom2, pos2, clip_or_wr)

                if index > 0:
                    event_pos, chrom, pos2, chrom2 = connect_left(all_aligns[index], all_aligns[index -1])
                    cigar_index = 0
                    # next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions, chrom2, pos2, pos2+1)
                    success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                           tell, cigar_index, event_pos, chrom2, pos2, clip_or_wr)

        if add_primary_link: # and paired_end:

            # Use mate information to generate cluster
            if clip_or_wr >= 2:  # within read  ==2
                chrom2 = chrom
                if r.cigartuples[cigar_index][0] != 1:  # not insertion, use length of cigar event
                    pos2 = cigar_pos2 #event_pos + r.cigartuples[cigar_index][1]
                else:
                    pos2 = event_pos

            if current_overlaps_roi and io_funcs.intersecter_int_chrom(overlap_regions, chrom2, pos2, pnext+1):
                # Probably too many reads in ROI to reliably separate out non-soft-clipped reads
                return

            if flag & 2 and not flag & 2048:  # Skip non-discordant/non-split that are same loci
                if abs(r.tlen) < clustering_dist and clip_or_wr < 2:  # != 2
                    return

            success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                   tell, cigar_index, event_pos, chrom2, pos2, clip_or_wr)

    else:  # Single end

        if clip_or_wr == 1:

            # Use SA tag to get chrom2 and pos2
            if r.has_tag("SA"):

                all_aligns, index = alignments_from_sa_tag(r, gettid)
                event = all_aligns[index]

                if len(all_aligns) == 1:
                    return  # shouldnt happen

                if index < len(all_aligns) - 1:  # connect to next
                    event_pos, chrom, pos2, chrom2 = connect_right(all_aligns[index], all_aligns[index + 1])
                    cigar_index = 0  #len(r.cigartuples) - 1

                    success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                   tell, cigar_index, event_pos, chrom2, pos2, clip_or_wr)

                if index > 0:
                    event_pos, chrom, pos2, chrom2 = connect_left(all_aligns[index], all_aligns[index -1])
                    cigar_index = len(r.cigartuples) - 1

                    success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                           tell, cigar_index, event_pos, chrom2, pos2, clip_or_wr)

        elif clip_or_wr >= 2:  # Sv within read
            chrom2 = r.rname
            if r.cigartuples[cigar_index][0] != 1:  # If not insertion
                pos2 = cigar_pos2  # event_pos + r.cigartuples[cigar_index][1]
            else:
                pos2 = event_pos

            success = add_to_graph(G, r, pe_scope, template_edges, node_to_name, genome_scanner, flag, chrom,
                                   tell, cigar_index, event_pos, chrom2, pos2, clip_or_wr)


cdef struct CigarEvent:
    int opp
    int cigar_index
    int event_pos
    int clip_or_wr
    int pos2
    int length
    bint cigar_skip


cdef CigarEvent make_cigar_event(int opp, int cigar_index, int event_pos, int clip_or_wr, int pos2, int length) nogil:
    cdef CigarEvent item
    item.opp = opp
    item.cigar_index = cigar_index
    item.event_pos = event_pos
    item.clip_or_wr = clip_or_wr
    item.pos2 = pos2
    item.length = length
    item.cigar_skip = 0
    return item


cpdef tuple construct_graph(genome_scanner, infile, int max_dist, int clustering_dist, int k=16, int m=7, int clip_l=21,
                            int min_sv_size=30,
                            int minimizer_support_thresh=2, int minimizer_breadth=3,
                            int minimizer_dist=10, int mapq_thresh=1, debug=None, int min_support=3, procs=1,
                            int paired_end=1):

    t0 = time.time()
    click.echo("Building graph with clustering distance {} bp, scope length {} bp".format(max_dist, clustering_dist),
               err=True)

    cdef TemplateEdges_t template_edges = TemplateEdges()  # Edges are added between alignments from same template, after building main graph
    cdef int event_pos, cigar_index, opp, length, added
    node_to_name = NodeToName()  # Map of nodes -> read ids


    # Infers long-range connections, outside local scope using pe information
    cdef PairedEndScoper_t pe_scope = PairedEndScoper(max_dist, clustering_dist, infile.header.nreferences)
    cdef long tell
    G = map_set_utils.Py_SimpleGraph()
    overlap_regions = genome_scanner.overlap_regions  # Get overlapper, intersect reads with intervals
    gettid = infile.gettid

    cdef int pos2, clip_or_wr #, n_aligned_bases
    # cdef bint del_was_last

    cdef vector[CigarEvent] events_to_add
    cdef vector[CigarEvent].iterator itr_events
    cdef CigarEvent v

    for chunk in genome_scanner.iter_genome():
        for r, tell in chunk:
            if r.mapq < mapq_thresh:
                continue

            pos2 = -1
            event_pos = r.pos
            added = 0

            # n_aligned_bases = 0
            # target_bases = 0
            # del_was_last = 0

            events_to_add.clear()
            if len(r.cigartuples) > 1:

                if r.has_tag("SA"):

                    # Set cigar-index to -1 means it is unset, will be determined during SA parsing
                    cigar_index = -1
                    clip_or_wr = 1  # Soft-clip or SA event
                    process_alignment(G, r, clip_l, max_dist, gettid,
                       overlap_regions, clustering_dist, pe_scope,
                       cigar_index, event_pos, paired_end, tell, genome_scanner,
                       template_edges, node_to_name,
                                      clip_or_wr, pos2)
                    added += 1

                for cigar_index, (opp, length) in enumerate(r.cigartuples):

                    if opp == 1:
                        if length >= min_sv_size:
                            clip_or_wr = 3  # Insertion type
                            pos2 = event_pos + length
                            events_to_add.push_back(make_cigar_event(opp, cigar_index, event_pos, clip_or_wr, pos2, length))
                            added += 1

                            # del_was_last = 0

                    elif opp == 2:
                        # Legacy code. Keep until nanopore has been tested
                        # Check if last deletion should be extended, short alignment block are merged with the del
                        # this only occurs until dist_to_last_sv has been reached
                        # if del_was_last and n_aligned_bases <= target_bases:
                        #
                        #     if length >= min_sv_size:
                        #         target_bases = min(150, max(100, int(length / 2)))
                        #
                        #         n_aligned_bases = 0
                        #         events_to_add.back().cigar_skip = 1
                        #
                        #     pos2 = event_pos + length
                        #     events_to_add.back().pos2 = pos2

                        if length >= min_sv_size:  # elif!
                            clip_or_wr = 2
                            pos2 = event_pos + length
                            events_to_add.push_back(make_cigar_event(opp, cigar_index, event_pos, clip_or_wr, pos2, length))

                            # target_bases = min(150, max(100, int(length / 2)))
                            # n_aligned_bases = 0

                            added += 1
                            # del_was_last = 1

                        event_pos += length

                    else:
                        if opp != 4 and opp != 5:
                            # if opp == 0:
                            #     n_aligned_bases += length
                            event_pos += length

            if not added:
                # Whole alignment will be used if mate information is available
                cigar_index = -1
                clip_or_wr = 1
                pos2 = -1
                process_alignment(G, r, clip_l, max_dist, gettid,
                           overlap_regions, clustering_dist, pe_scope,
                           cigar_index, event_pos, paired_end, tell, genome_scanner,
                           template_edges, node_to_name,
                                  clip_or_wr, pos2)

            if not events_to_add.empty():
                itr_events = events_to_add.begin()

                while itr_events != events_to_add.end():
                    v = dereference(itr_events)
                    if v.cigar_skip:
                        pos2 = v.pos2
                    else:
                        pos2 = v.event_pos + v.length  # fall back on original cigar event length
                    # echo(f"{v.event_pos}-{v.pos2}", v.cigar_index)
                    process_alignment(G, r, clip_l, max_dist, gettid,
                               overlap_regions, clustering_dist, pe_scope,
                               v.cigar_index, v.event_pos, paired_end, tell, genome_scanner,
                               template_edges, node_to_name,
                                      v.clip_or_wr, v.pos2)

                    preincrement(itr_events)
            # echo("---")

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
    #
    # echo("parts", partitions)
    # echo("s_between", sb)
    # echo("s_within", support_within)
    # quit()
    # echo("n2n", n2n.keys())

    return {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}

