#cython: language_level=3

import click
from collections import deque, defaultdict

from libcpp.vector cimport vector as cpp_vector
from libcpp.deque cimport deque as cpp_deque
from libcpp.pair cimport pair as cpp_pair
from libcpp.string cimport string as cpp_string

import cython

from libc.stdlib cimport abs as c_abs
from libc.math cimport fabs as c_fabs

# from pysam.libcalignmentfile cimport AlignmentFile
# from pysam.libcalignedsegment cimport AlignedSegment
# from pysam.libchtslib cimport bam1_t, BAM_CIGAR_SHIFT, BAM_CIGAR_MASK
from libc.stdint cimport uint32_t

from dysgu.map_set_utils cimport hash as xxhasher

# ctypedef cpp_vector[int] int_vec_t
# ctypedef cpp_pair[int, int] get_val_result

# ctypedef Py_Int2IntVecMap[int, int_vec_t] node_dict_t
# ctypedef Py_IntVec2IntMap[int_vec_t, int] node_dict2_r_t

import mmh3

ctypedef cpp_pair[int, int] cpp_item
ctypedef cpp_pair[long, int] cpp_long_item


def echo(*args):
    click.echo(args, err=True)


cdef class Py_DiGraph:
    """DiGraph, weighted"""
    def __cinit__(self):
        self.thisptr = new DiGraph()
    def __dealloc__(self):
        del self.thisptr
    cdef int addNode(self):
        return self.thisptr.addNode()
    cdef int hasEdge(self, int u, int v):
        return self.thisptr.hasEdge(u, v)
    cdef void addEdge(self, int u, int v, int w):
        self.thisptr.addEdge(u, v, w)
    cdef void updateEdge(self, int u, int v, int w):
        self.thisptr.addEdge(u, v, w)
    cdef int numberOfNodes(self) nogil:
        return self.thisptr.numberOfNodes()
    cdef cpp_vector[cpp_pair[int, int]] forInEdgesOf(self, int u) nogil:
        return self.thisptr.forInEdgesOf(u)
    cdef cpp_vector[int] neighbors(self, int u) nogil:
        return self.thisptr.neighbors(u)
    cdef float node_path_quality(self, int u, int v, int w)  nogil:
        return self.thisptr.node_path_quality(u, v, w)


cdef class Py_SimpleGraph:
    """Graph, weighted"""
    def __cinit__(self):
        self.thisptr = new SimpleGraph()
    def __dealloc__(self):
        del self.thisptr
    cpdef int addNode(self):
        return self.thisptr.addNode()
    cpdef int hasEdge(self, int u, int v):
        return self.thisptr.hasEdge(u, v)
    cpdef void addEdge(self, int u, int v, int w):
        self.thisptr.addEdge(u, v, w)
    cpdef int edgeCount(self):
        return self.thisptr.edgeCount()
    cpdef int weight(self, int u, int v):
        return self.thisptr.weight(u, v)
    cpdef cpp_vector[int] neighbors(self, int u):
        return self.thisptr.neighbors(u)
    cpdef void removeNode(self, int u):
        self.thisptr.removeNode(u)
    cpdef cpp_vector[int] connectedComponents(self):
        return self.thisptr.connectedComponents()
    cpdef int showSize(self):
        return self.thisptr.showSize()


cdef class Py_Int2IntMap:
    """Fast integer to integer unordered map using robin_hood flat map"""
    def __cinit__(self):
        self.thisptr = new Int2IntMap()
    def __dealloc__(self):
        del self.thisptr
    cdef void insert(self, int key, int value) nogil:
        self.thisptr.insert(key, value)
    cdef void erase(self, int key) nogil:
        self.thisptr.erase(key)
    cdef int has_key(self, int key) nogil:
        return self.thisptr.has_key(key)
    cdef int get(self, int key) nogil:
        return self.thisptr.get(key)
    cdef get_val_result get_value(self, int key) nogil:
        return self.thisptr.get_value(key)
    cdef int size(self) nogil:
        return self.thisptr.size()


cdef class Py_IntSet:
    """Fast set using robin_hood unordered set"""
    def __cinit__(self):
        self.thisptr = new IntSet()
    def __dealloc__(self):
        del self.thisptr
    cdef void insert(self, int key) nogil:
        self.thisptr.insert(key)
    cdef void erase(self, int key) nogil:
        self.thisptr.erase(key)
    cdef int has_key(self, int key) nogil:
        return self.thisptr.has_key(key)
    cdef int size(self) nogil:
        return self.thisptr.size()


cdef int cigar_exists(r):
    if r.cigartuples:
        return 1
    return 0


cdef tuple clip_sizes(r):
    c = r.cigartuples
    if not c:
        return 0, 0

    cdef int left = 0
    cdef int right = 0

    if c[0][0] == 4:
        left = c[0][1]
    if c[-1][0] == 4:
        right = c[-1][1]
    return left, right


cdef tuple clip_sizes_hard(r):
    c = r.cigartuples
    if not c:
        return 0, 0

    cdef int left = 0
    cdef int right = 0
    c1 = c[0][0]
    if c1 == 4 or c1 == 5:
        left = c[0][1]
    c1 = c[-1][0]
    if c1 == 4 or c1 == 5:
        right = c[-1][1]
    return left, right


cdef int cigar_clip(r, int clip_length):

    c = r.cigartuples
    if not c:
        return 0
    if (c[0][0] == 4 and c[0][1] >= clip_length) or (c[-1][0] == 4 and c[-1][1] >= clip_length):
        return 1
    return 0


cdef int is_overlapping(int x1, int x2, int y1, int y2) nogil:
    return int(max(x1, y1) <= min(x2, y2))


cdef bint is_reciprocal_overlapping(int x1, int x2, int y1, int y2) nogil:
    # Insertions have same x1/y1 position, use another measure
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

    position_distance = c_fabs(center1 - center2) # 1 #distance_normalizer
    if position_distance > 2000:
        return 0
    max_span = max(span1, span2)
    span_distance = <float>c_abs(span1 - span2) / max_span
    # echo("pd", position_distance, center1, center2)
    # echo((position_distance / max_span), span_distance, center1, center2 )
    if (position_distance / max_span) < 0.2 and span_distance < 0.3:
    # if position_distance < 100 and span_distance < 0.08:
        return 1
    return 0


cdef float position_distance(int x1, int x2, int y1, int y2) nogil:
    # https://github.com/eldariont/svim/blob/master/src/svim/SVIM_clustering.py
    cdef int span1, span2
    cdef float center1, center2
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
    return c_fabs(center1 - center2)


cdef void sliding_window_minimum(int k, int m, str s, unordered_set[long]& found):
    """End minimizer. A iterator which takes the size of the window, `k`, and an iterable,
    `li`. Then returns an iterator such that the ith element yielded is equal
    to min(list(li)[max(i - k + 1, 0):i+1]).
    Each yield takes amortized O(1) time, and overall the generator takes O(k)
    space.
    https://github.com/keegancsmith/Sliding-Window-Minimum/blob/master/sliding_window_minimum.py"""

    # Cpp version
    cdef int i = 0
    cdef int end = len(s) - m + 1
    cdef int last_idx = end - 1
    cdef cpp_deque[cpp_long_item] window2
    cdef long int hx2

    cdef bytes s_bytes = bytes(s.encode("ascii"))

    # cdef char* my_ptr #= <char*>&my_view[0]
    # cdef const unsigned char[:] sub
    cdef const unsigned char* sub_ptr = s_bytes
    # cdef cpp_set[int] seen2  # Using
    # cdef cpp_u_set[int] seen2
    # seen2 = set([])

    # cdef cpp_item last

    with nogil:

        for i in range(end):
            # xxhasher(bam_get_qname(r._delegate), len(qname), 42)

            # kmer = s[i:i+m]
            # if kmer == len(s) * kmer[0]:
            #     continue  # skip homopolymers

            # hx2 = mmh3.hash(kmer, 42)
            hx2 = xxhasher(sub_ptr, m, 42)
            if i == 0 or i == last_idx:
                found.insert(hx2)
                # seen2.add(hx2)
            # elif i == end - 1:
                # seen2.add(hx2)
                # break

            # sub = s_bytes[i:i+m]
            # sub = b'TGGAGAAGAG'


            # sub_ptr = sub
            # sub = view[i:i+m]

            # hx2 = xxhasher(sub_ptr, len(s_bytes), 42)
            # echo(sub_ptr[0], hx2)
            # echo(mmh3.hash(s[i:i+m], 42), s_bytes[i:i+m], hx2)
            # quit()
            while window2.size() != 0 and window2.back().first >= hx2:
                window2.pop_back()

            window2.push_back(cpp_long_item(hx2, i))
            while window2.front().second <= i - k:
                window2.pop_front()

            # i += 1

            sub_ptr += 1

            minimizer_i = window2.front().first

            #if minimizer_i not in seen2:
            # seen2.add(minimizer_i)

            found.insert(minimizer_i)
            # if seen2.find(minimizer_i) == seen2.end():
            #     seen2.insert(minimizer_i)

    # return seen2  #set(seen2)


cdef str left_soft_clips(str seq, int code_length):
    return seq[0:code_length]


cdef str right_soft_clips(str seq, int code_length):
    return seq[len(seq) - code_length:]


class ClipScoper:
    """Keeps track of which reads are in scope. Maximum distance depends on the template insert_median"""
    def __init__(self, int max_dist, int k, int m, int clip_length, int minimizer_support_thresh,
                 int minimizer_breadth, int read_length):
        self.max_dist = max_dist
        self.k = k
        self.w = m
        self.clip_length = clip_length

        self.scope_left = deque([])
        self.scope_right = deque([])
        self.clip_table = {0: defaultdict(set), 1: defaultdict(set)}
        self.read_minimizers = defaultdict(set)

        self.current_chrom = 0
        self.minimizer_support_thresh = minimizer_support_thresh
        self.minimizer_matches = minimizer_breadth

        self.target_density = 2. / (m + 1)
        self.upper_bound_n_minimizers = read_length * self.target_density

    def _add_m_find_candidates(self, clip_seq, int name, int idx, int position):

        cdef unordered_set[long] clip_minimizers
        sliding_window_minimum(self.k, self.w, clip_seq, clip_minimizers)
        if clip_minimizers.empty():
            return

        # add read minimizers from read, and find partners
        # idx 0 is for left clips, 1 is for right clips
        target_counts = defaultdict(int)
        total_m_found = 0
        cdef int find_candidate = 1
        cdef long int m
        cdef long int item
        res = set([])

        n_local_minimizers = len(self.clip_table[idx])
        n_local_reads = len(self.scope_left) if idx == 0 else len(self.scope_right)

        upper_bound = (1 + (n_local_reads * 0.05)) * self.upper_bound_n_minimizers
        # echo(n_local_minimizers, self.upper_bound_n_minimizers, upper_bound)
        if n_local_minimizers > upper_bound:
            # echo(n_local_reads, n_local_minimizers, upper_bound)
            find_candidate = 0

        for m in clip_minimizers:

            min_id = (m, idx)
            self.read_minimizers[name].add(min_id)
            minimize_table = self.clip_table[idx]

            if m not in minimize_table:
                minimize_table[m].add((position, name))
                continue

            elif find_candidate:  # Look for suitable partners
                targets = minimize_table[m]
                for item_position, item in targets:
                    if abs(item_position - position) < 10:

                        target_counts[item] += 1
                        total_m_found += 1
                        support = (total_m_found / 2) + len(target_counts)
                        # echo(support)
                        if support >= self.minimizer_support_thresh:
                            res.update(target_counts.keys())

                        if len(res) >= 4:  # Maximum edges for each read
                            find_candidate = 0
                            break

            minimize_table[m].add((position, name))

        return res

    def _refresh_scope(self, scope, position):
        # Remove out of scope reads and minimizers
        # left
        while True:
            if len(scope) == 0:
                break
            if abs(scope[0][0] - position) > self.max_dist:
                name = scope.popleft()  # position, input_read
                if name[1] in self.read_minimizers:
                    for m, idx in self.read_minimizers[name[1]]:
                        minimizer_table = self.clip_table[idx]
                        if m in minimizer_table:
                            minimizer_table[m].remove(name)
                            if len(minimizer_table[m]) == 0:
                                del minimizer_table[m]
                    del self.read_minimizers[name[1]]
            else:
                break

    def _insert(self, str seq, int cigar_start, int cigar_end, int input_read, int position):
        # Find soft-clips of interest
        cdef set clip_set
        cdef str clip_seq
        targets = set([])
        if cigar_start >= self.clip_length:
            self._refresh_scope(self.scope_left, position)

            clip_seq = left_soft_clips(seq, cigar_start) #"".join(left_soft_clips(seq, cigar_start)[::-1])
            targets |= self._add_m_find_candidates(clip_seq, input_read, 0, position)

            self.scope_left.append((position, input_read))

        if cigar_end >= self.clip_length:
            self._refresh_scope(self.scope_right, position)

            clip_seq = right_soft_clips(seq, cigar_end)
            # clip_set = sliding_window_minimum(self.k, self.w, clip_seq)
            targets |= self._add_m_find_candidates(clip_seq, input_read, 1, position)

            self.scope_right.append((position, input_read))

        return targets

    def update(self, int input_read, str seq, int cigar_start, int cigar_end, int chrom, int position):

        # if len(self.scope) == 0:
        #     self.scope.append((position, input_read))
        #     self.current_chrom = chrom
        #     self._insert(seq, cigar_start, cigar_end, input_read, position)  # Add minimizers
        #     return
        #
        # elif chrom != self.current_chrom:
        #     # Empty scope on new chromosome
        #     self.scope = deque([(position, input_read)])
        #     # self.minimizer_table = defaultdict(set)
        #     self.clip_table = {0: defaultdict(set), 1: defaultdict(set)}
        #
        #     self.read_minimizers = defaultdict(set)
        #     self.current_chrom = chrom
        #     self._insert(seq, cigar_start, cigar_end, input_read, position)
        #     return

        if chrom != self.current_chrom:
            # Empty scope on new chromosome
            self.scope_left = deque([])
            self.scope_right = deque([])
            self.clip_table = {0: defaultdict(set), 1: defaultdict(set)}
            self.read_minimizers = defaultdict(set)
            self.current_chrom = chrom

        d = self._insert(seq, cigar_start, cigar_end, input_read, position)
        if len(d) > 0:
            return d  # Best candidate with most overlapping minimizers
