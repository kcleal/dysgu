#cython: language_level=3, boundscheck=False

from __future__ import absolute_import
import time
import numpy as np
import click
from collections import defaultdict, deque
import itertools
import mmh3
import sortedcontainers
import operator
import ncls
import networkit as nk
import pysam
cimport numpy as c_np

from libcpp.vector cimport vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
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

nk.setSeed(100, True)

ctypedef cpp_pair[int, int] cpp_item
ctypedef map_set_utils.Py_IntSet Py_IntSet


def echo(*args):
    click.echo(args, err=True)


cdef class Table:
    cdef vector[c_np.int64_t] starts
    cdef vector[c_np.int64_t] ends
    cdef vector[c_np.int64_t] values

    cpdef void add(self, int s, int e, int v):
        with nogil:
            self.starts.push_back(s)
            self.ends.push_back(e)
            self.values.push_back(v)

    def get_val(self, v):
        cdef vector[c_np.int64_t] values = v
        cdef c_np.ndarray[c_np.int64_t] a = np.empty(values.size(), dtype=np.int)
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
    cdef cpp_deque[cpp_item] window2
    cdef int hx2
    # cdef cpp_set[int] seen2  # Using
    # cdef cpp_u_set[int] seen2
    seen2 = set([])

    # cdef cpp_item last
    for i in range(0, len(s) - m + 1):

        hx2 = int(mmh3.hash(s[i:i+m], 42))

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


# cdef set sliding_window_minimum_p(int k, int m, str s):
#     '''
#     End minimizer. A iterator which takes the size of the window, `k`, and an iterable,
#     `li`. Then returns an iterator such that the ith element yielded is equal
#     to min(list(li)[max(i - k + 1, 0):i+1]).
#     Each yield takes amortized O(1) time, and overall the generator takes O(k)
#     space.
#     https://github.com/keegancsmith/Sliding-Window-Minimum/blob/master/sliding_window_minimum.py
#     '''
#
#     # Python version for posterity
#     window = deque()
#     cdef int i = 0
#     cdef set seen = set([])
#     cdef int hx
#     for i in range(0, len(s) - m + 1):
#
#         hx = mmh3.hash(s[i:i+m], 42)
#
#         while window and window[-1][0] >= hx:
#             window.pop()
#
#         window.append((hx, i))
#
#         while window[0][1] <= i - k:
#             window.popleft()
#
#         i += 1
#
#         minimizer_i = window[0][0]
#         if minimizer_i not in seen:
#             seen.add(minimizer_i)
#     return seen


class ClipScoper(object):
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

    def _add_m_find_candidates(self, set clip_set, int name, int idx):
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

            elif find_candidate:
                targets = self.minimizer_table[min_id]

                for item in targets:

                    target_counts[item] += 1
                    # if name == 63018:
                    #     echo(item, target_counts)
                    # if target_counts[item] >= self.minimizer_matches and \
                    #         len(target_counts) >= self.minimizer_support_thresh:
                    #     res.add(item)
                    # if len(res) >= 3:
                    #     find_candidate = 0
                    #     break
                    if target_counts[item] >= self.minimizer_matches and \
                            len(target_counts) >= self.minimizer_support_thresh:
                        res.add(item)
                    if len(res) >= 3:  # Maximum of 3 edges for each read
                        find_candidate = 0
                        break

                            # Do some testing to see if seqs match
                        # find_candidate = 0

            self.minimizer_table[min_id].add(name)
        # if name in {46411}:
        #     echo("Res", name, res, target_counts)
            # echo(self.minimizer_table)
            # echo(self.read_minimizers)
            # echo()
            # quit()

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
            # if input_read in {63010, 63018}:
            #     echo(input_read, sorted(clip_set))
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
        # debug = {8, 429, 475, 64, 7}
        # if input_read in debug:
        #     echo("Debugging", input_read, position, "---------")
        #     echo(len(self.scope))
        #     echo("Scope", list(self.scope)[:100], input_read, position)
        #     echo(self.max_dist, (self.current_chrom, chrom, position))
        #     echo()
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
        # if input_read in debug:
        #     echo("Scope2", list(self.scope)[:100])

        cdef set d = self._insert(seq, cigar_start, cigar_end, input_read)
        # if input_read in debug:
        #     echo("Set", d)
            # quit()
        if len(d) > 0:
            return d  # max(d, key=d.get)  # Best candidate with most overlapping minimizers


class ForwardScope:

    def __init__(self):
        self.inscope = 0
        self.data = sortedcontainers.SortedKeyList(key=operator.itemgetter(0))
    def add(self, v):
        self.data.add(v)
        self.inscope += 1
    def drop_one(self):
        # Keep a reference count of number in local scope, means multiple values can be deleted by reinitializing list
        self.inscope -= 1
        if self.inscope == 0:
            self.data = sortedcontainers.SortedKeyList(key=operator.itemgetter(0))
        elif self.inscope < 0:
            raise ValueError
    def empty(self):
        self.inscope = 0
        self.data = sortedcontainers.SortedKeyList(key=operator.itemgetter(0))


class PairedEndScoper:

    def __init__(self, max_dist):
        self.clst_dist = 3*max_dist  # Keep reads within this distance in scope  # todo sensible threshold
        self.max_dist = 3*max_dist
        self.local_chrom = None
        self.local_scope = deque([])
        self.forward_scope = defaultdict(lambda: ForwardScope())


    def update(self, node_name, current_chrom, current_pos, chrom2, pos2):

        if current_chrom != self.local_chrom:

            self.local_chrom = current_chrom
            self.local_scope = deque([(current_pos, chrom2)])

            self.forward_scope = defaultdict(lambda: ForwardScope())
            self.forward_scope[chrom2].add((pos2, node_name, current_pos))

            return None

        while len(self.local_scope) > 0:
            if abs(self.local_scope[0][0] - current_pos) > self.clst_dist:
                local_p, f_chrom = self.local_scope.popleft()
                # echo("here")
                self.forward_scope[f_chrom].drop_one()

            else:
                break

        # Find other nodes in forward_scope
        t = self.forward_scope[chrom2].data
        len_t = len(t)
        if len_t == 0:
            t.add((pos2, node_name, current_pos))
            return None

        # pnext - max_dist, pnext + max_dist
        idx = t.bisect_left((pos2, None))
        found = []
        if idx == 0:
            # echo("t", t, idx)
            item_forward_pos, item_node, item_local_pos = t[idx]

            if abs(item_local_pos - current_pos) < self.max_dist and \
               abs(item_forward_pos - pos2) < self.max_dist:
                    found.append(item_node)

        if idx == len_t:
            item_forward_pos, item_node, item_local_pos = t[idx - 1]

            if abs(item_local_pos - current_pos) < self.max_dist and \
               abs(item_forward_pos - pos2) < self.max_dist:
                    found.append(item_node)

        else:  # Add upstream and downstream
            item_forward_pos, item_node, item_local_pos = t[idx]

            if abs(item_local_pos - current_pos) < self.max_dist and \
               abs(item_forward_pos - pos2) < self.max_dist:
                    found.append(item_node)

            item_forward_pos, item_node, item_local_pos = t[idx - 1]

            if abs(item_local_pos - current_pos) < self.max_dist and \
               abs(item_forward_pos - pos2) < self.max_dist:
                    found.append(item_node)

        # Add to scopes
        self.local_scope.append((current_pos, chrom2))
        self.forward_scope[chrom2].add((pos2, node_name, current_pos))
        return found



def construct_graph(genome_scanner, infile, int max_dist, int clustering_dist, int k=16, int m=7, int clip_l=21,
                    int minimizer_support_thresh=2, int minimizer_breadth=3,
                    int minimizer_dist=10, debug=None, int min_support=3, procs=1):

    t0 = time.time()
    click.echo("Building graph with clustering distance {}bp".format(max_dist), err=True)

    all_flags = defaultdict(list)  # Linking clusters together rname: set([(flag, pos)..])

    node_to_name = []

    scope = ClipScoper(minimizer_dist, k=k, m=m, clip_length=clip_l,
                       minimizer_support_thresh=minimizer_support_thresh,
                       minimizer_breadth=minimizer_breadth)

    pe_scope = PairedEndScoper(max_dist)

    overlap_regions = genome_scanner.overlap_regions

    cdef int count = 0
    cdef int buf_del_index = 0
    cdef int ii = 0
    cdef int grey_added = 0

    cdef int flag, pos, rname, pnext, rnext, node_name, node_name2, chrom, chrom2, clip_left, clip_right
    cdef str qname, seq

    cdef int warn = 0
    cdef int ol_include, add_primary_link
    cdef int current_overlaps_roi, next_overlaps_roi

    t0 = time.time()

    G = nk.graph.Graph(weighted=True)

    # testout = pysam.AlignmentFile("/Users/kezcleal/Documents/Data/fusion_finder_development/AshkenazimTrio/calls/strange.bam", "wb", template=infile)
    # cdef list chunk
    # cdef AlignedSegment r

    debug = "D00360:18:H8VC6ADXX:2:2210:15612:70818"
    debug_nodes = set([])
    for chunk in genome_scanner.iter_genome():
        for r, tell in chunk:

            # testout.write(r)

            seq = r.seq

            clip_left, clip_right = map_set_utils.clip_sizes(r)

            chrom = r.rname
            pos = r.pos
            flag = r.flag
            qname = r.qname

            n1 = (mmh3.hash(qname, 42), flag, pos, chrom, tell)
            node_name = G.addNode()
            node_to_name.append(n1)

                # quit()
            genome_scanner.add_to_buffer(r, node_name)

            all_flags[qname].append(node_name)  # Used for paired-end edges

            if qname == debug:
                echo(qname, n1, "NODE NAME", node_name, clip_left, clip_right)
                debug_nodes.add(node_name)

            # Cluster soft-clips
            # if cigartuples[0][0] == 4 or cigartuples[-1][0] == 4:
            if clip_left > clip_l or clip_right > clip_l:
                best_candidates = scope.update(node_name, seq, clip_left, clip_right, chrom, pos)

                # if qname == debug:
                #     echo(qname, n1, "NODE NAME", node_name)
                #     echo(best_candidates)
                #     debug_nodes.add(node_name)

                if best_candidates is not None:
                    for node_name2 in best_candidates:
                        if not G.hasEdge(node_name, node_name2):
                            G.addEdge(node_name, node_name2, w=3)  # weight 3 for black edge
                            # if node_name in debug_nodes or node_name2 in debug_nodes:
                            #     echo("black edge", node_name, node_name2)
                            #     echo(node_name, str(r.qname))

            #if flag & 2 and not flag & 2048:  # Skip non-discordant. Supplementary have rnext:pnext mixed up
            #    if abs(r.tlen) < clustering_dist:
            #        continue

            # Skip include regions, coverage too high to work with
            # if tree is not None and data_io.intersecter(tree, infile.get_reference_name(r.rname), pos, pos + 1):
            #     continue

            # Cluster paired-end mates
            pnext = r.pnext
            rnext = r.rnext

            # Special treatment of supplementary and local reads; need to decide where the partner is
            # Either use the rnext:pnext or a site(s) listed in the SA tag: The rnext:pext can often be at the
            # same loci as the read which leads to problems when linking no.2 edges
            add_primark_link = 1

            current_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions,
                                                                    r.rname,
                                                                    pos,
                                                                    pos+1)
            # if node_name in debug_nodes:
            #     echo("current overlaps roi", node_name, current_overlaps_roi)

            if chrom == rnext and abs(pnext - pos) < (max_dist * 3):  # Same loci

                if r.has_tag("SA"):  # Parse SA, first alignment is the other read primary line
                    sa = r.get_tag("SA").split(",", 2)
                    chrom2 = infile.gettid(sa[0])
                    pos2 = int(sa[1])
                    if chrom2 != chrom or abs(pos2 - pos) > (max_dist * 3):

                        add_primark_link = 0
                        next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions,
                                                                             chrom2,
                                                                             pos2,
                                                                             pos2+1)

                        # if not next_overlaps_roi:
                        if current_overlaps_roi and next_overlaps_roi:
                            continue
                        else:
                            other_nodes = pe_scope.update(node_name, chrom, pos, chrom2, pos2)

                            # if qname == debug:
                            #     echo(node_name, other_nodes)

                            if other_nodes:
                                for other_node in other_nodes:

                                    if not G.hasEdge(node_name, other_node):
                                        G.addEdge(node_name, other_node, w=2)

                                        # if qname == debug:
                                        #     echo(node_name, other_node, "SA")

            #
            if add_primark_link == 1:
                next_overlaps_roi = io_funcs.intersecter_int_chrom(overlap_regions,
                                                                     rnext,
                                                                     pnext,
                                                                     pnext+1)

                if current_overlaps_roi and next_overlaps_roi:
                    continue

                if flag & 2 and not flag & 2048:  # Skip non-discordant that are same loci
                    if abs(r.tlen) < clustering_dist:
                       continue

                other_nodes = pe_scope.update(node_name, chrom, pos, rnext, pnext)
                if other_nodes:
                    for other_node in other_nodes:
                        # pass
                        if not G.hasEdge(node_name, other_node):
                            G.addEdge(node_name, other_node, w=2)

                            # if node_name in debug_nodes or node_name2 in debug_nodes:
                            #     echo("normal edge", node_name, node_name2, current_overlaps_roi, next_overlaps_roi)
                            #     echo(node_name, str(r.qname))

    # testout.close()
    # import subprocess
    # subprocess.call("samtools index /Users/kezcleal/Documents/Data/fusion_finder_development/AshkenazimTrio/calls/strange.bam", shell=True)

    click.echo(f"Processed minimizers {time.time() - t0}", err=True)
    # Make nested containment list

    # t1 = time.time()
    # nc = {rnext: val.containment_list() for rnext, val in interval_table.items()}
    # click.echo(f"Containement list made in {time.time() - t1}", err=True)
    # click.echo(f"Length of node_to_next {len(node_to_next)}", err=True)
    # t1_5 = time.time()
    # cdef int u, v
    # cdef set grey_edges = set([])
    # cdef int start, end
    #
    # for u, chrom, pos, rnext, pnext in node_to_next:
    #     count = 0
    #
    #     for start, end, v in (nc[rnext].find_overlap(pnext - max_dist, pnext + max_dist)):
    #         if u != v and not G.hasEdge(u, v):
    #
    #             other_node = node_to_name[v]
    #             other_pos = other_node[2]
    #             other_chrom = other_node[3]
    #
    #             if other_chrom == chrom and abs(other_pos - pos) < (2*max_dist):
    #                 G.addEdge(u, v, w=2)
    #                 count += 1
    #                 if count > 4:
    #                     break
    #
    # click.echo(f"Containement list run in {time.time() - t1_5}", err=True)

    # echo(">", all_flags[debug])

    t2 = time.time()
    cdef list node_names
    for qname, node_names in all_flags.items():  # 2 reads, or 3 if supplementary read
        for u, v in itertools.combinations_with_replacement(node_names, 2):  # Check all edges within template
            if u != v and not G.hasEdge(u, v):
                G.addEdge(u, v, w=1)

    click.echo(f"Paired end info in {time.time() - t2}", err=True)

    cdef dict component

    read_buffer = genome_scanner.read_buffer

    click.echo("Constructed graph, processing block model", err=True)
    for component in get_block_components(G, node_to_name, infile, max_dist, minimizer_dist, read_buffer, min_support,
                                          debug_nodes):

        yield component


def get_block_components(G, list node_to_name, infile, int max_dist, int read_length, dict read_buffer,
                         int min_support, debug_nodes):

    # Turn graph into a block model, first split into connected components,
    # Then for each component, split into block nodes which are locally interacting nodes (black and grey edges)
    # block nodes edges result from collapsing white edges
    cc = nk.components.ConnectedComponents(G)
    cc.run()

    cdef dict partitions, support_within, reads
    cdef int v

    big_components = []
    t0 = time.time()

    cdef int cnt = 0

    for component in cc.getComponents():

        cnt += 1

        if len(component) >= min_support:
            res = proc_component(node_to_name, component, read_buffer, infile, G, min_support)
            # if any(itemn in debug_nodes for itemn in component):
            #     echo("hi len component", len(component))
            #     echo(res)
            #     for qname, align in res["reads"].items():
            #         echo(qname, str(align).replace("\t", " "))

            if res:
                yield res

        # Reduce size of graph
        for v in component:
            G.removeNode(v)

    click.echo(f"Processed components n={cnt} {time.time() - t0}", err=True)


cpdef dict get_reads(infile, dict sub_graph_reads):

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


cdef dict proc_component(list node_to_name, list component, dict read_buffer, infile, G, int min_support):

    n2n = {}
    reads = {}

    cdef int v
    for v in component:
        if v in read_buffer:
            reads[v] = read_buffer[v]
        else:
            n2n[v] = node_to_name[v]


    # Get any reads that are in read buffer, leave the rest in n2n (these will be fetched later)


    # if len(missing) > 0:
    #     click.echo(f"Warning missing reads n={len(missing)}", err=True)
    #     echo(missing)
    #     for item in missing:
    #         del n2n[item]
    #     component = list(set(component).difference(missing))
    #     for rn in missing:
    #         G.removeNode(rn)
    #     if len(component) == 0:
    #         return None


    # Explore component for locally interacting nodes; create partitions using these
    partitions = get_partitions(G, component)

    support_between, support_within = count_support_between(G, partitions)

    # echo("support wthi", support_within)
    # for k, v in partitions.items():
    #     if 580 in v:
    #         echo("in parts")
    #         echo(partitions)
    #         echo(support_between)
    #         echo(support_within)
    #         quit()

    if len(support_between) == 0 and len(support_within) == 0:
        return {"parts": [], "s_between": [], "reads": reads, "s_within": [], "n2n": n2n}

    sb = {}
    kept = set([])
    deleted = set([])
    for edge, vd in support_between.items():

        sup = sum([len(vv) for vv in vd.values()])
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

    # debug_component(component, node_to_name, support_between, support_within, G, partitions, {157240, 157241, 157242})
    return {"parts": partitions, "s_between": sb, "reads": reads, "s_within": support_within, "n2n": n2n}


cdef list subgraph_from_nodes(G, list nodes):
    # Mark all the vertices as not visited
    edges_found = set([])
    cdef int u, v
    e = []
    for u in nodes:
        for v in G.neighbors(u):
            if (u, v) in edges_found or (v, u) in edges_found:
                continue
            e.append((u, v, {"w": G.weight(u, v)}))
            edges_found.add((u, v))
    return e


cdef tuple BFS_local(G, int source):

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
                if G.weight(u, v) > 1:

                    if u not in nodes_found:
                        nodes_found.add(u)
                    if v not in nodes_found:
                        nodes_found.add(v)
                        queue.append(v)

        visited.add(u)

    return nodes_found, visited


cdef dict get_partitions(G, list nodes):

    seen = set([])  # Dont use Py_IntSet, unions are not quick at the moment

    cdef int u, v, i
    parts = []
    for u in nodes:
        if u in seen:
            continue

        for v in G.neighbors(u):
            if v in seen:
                continue

            if G.weight(u, v) > 1:  # 2 or 3
                found, visited_local = BFS_local(G, u)
                seen |= visited_local
                if found:
                    parts.append(found)
        seen.add(u)
    return {i: p for i, p in enumerate(parts)}


cpdef tuple count_support_between(G, dict parts):
    # todo sometimes supplementary reads get confused with support within (use SA tag to get mate info)
    cdef int i, j, node, child, any_out_edges
    cdef set p
    cdef tuple t

    if len(parts) == 0:
        return {}, {}
    elif len(parts) == 1:
        return {}, {list(parts.keys())[0]: list(parts.values())[0]}

    # Make a table to count from
    p2i = dict()
    for i, p in parts.items():
        p2i.update({node: i for node in p})


    # Count the links between partitions. Split reads into sets ready for assembly
    # counts (part_a, part_b): {part_a: [node 1, node 2 ..], part_b: [node4, ..] }
    counts = {}

    # self_counts = dict()
    self_counts = defaultdict(list)

    # seen = set([])
    cdef Py_IntSet seen_nodes = map_set_utils.Py_IntSet()

    for i, p in parts.items():

        for node in p:
            any_out_edges = 0  # Keeps track of number of outgoing pairs, or self edges

            # if node in seen:
            #     continue
            if seen_nodes.has_key(node):
                continue

            for child in G.neighbors(node):

                if not child in p2i or seen_nodes.has_key(child):
                    continue  # Exterior child, not in any partition
                j = p2i[child]

                if j != i:
                    any_out_edges = 1
                    if j < i:
                        t = (j, i)
                    else:
                        t = (i, j)

                    if t not in counts:
                        counts[t] = defaultdict(set)

                    counts[t][i].add(node)
                    counts[t][j].add(child)

            # seen.add(node)
            seen_nodes.insert(node)
            # Count self links, important for resolving small SVs
            if any_out_edges == 0: #G.weight(node, child) == 2:

                # if i in self_counts and node in self_counts[i]:
                #     raise ValueError
                self_counts[i].append(node)

                # if i not in self_counts:
                #     self_counts[i] = 1
                # else:
                #     self_counts[i] += 1

    return counts, self_counts


def debug_component(component, node_to_name, support_between, support_within, G, partitions, targets):

    for cmp in component:
        if cmp in targets:
            echo("Supportbetween", support_between, "Support within", support_within)
            import networkx as nx
            nxG = nx.Graph()
            nxG.add_edges_from(subgraph_from_nodes(G, component))
            nx.write_gml(nxG, "/Users/kezcleal/Documents/Data/fusion_finder_development/NA12878/giab_compare/alignments/testnewtest.gml")
            echo("Partitons")
            for pp1, ppv in partitions.items():
                echo(pp1, ppv)
            echo()
            for cmp2 in sorted(component):
                echo(cmp2, node_to_name[cmp2])
            quit()
            break
