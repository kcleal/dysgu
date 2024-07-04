#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False
#distutils: language=c++
from libcpp.vector cimport vector
from libcpp.algorithm cimport lower_bound
from libc.stdlib cimport abs as c_abs


cdef bint is_overlapping(int x1, int x2, int y1, int y2):
    return max(x1, y1) <= min(x2, y2)


cdef struct Interval:
    int end
    int covered_end


cdef class IntervalSet:
    """An interval set for searching with sorted reference intervals and optionally sorted query points/intervals

    Reference intervals must be added in sorted order
    If query points are not sorted, a binary search is used, otherwise a linear search is used

    :param with_data: Whether to record additional data along with an interval
    :type with_data: bool
    :return: IntervalSet class for performing intersections
    :rtype: IntervalSet
    """
    cdef vector[int] starts
    cdef vector[Interval] ends
    cdef int last_r_start, last_q_start
    cdef int current_r_end, current_r_start
    cdef size_t index
    cdef bint add_data, bool_only
    cdef list data
    cdef int distance_threshold
    def __init__(self, with_data=False, bool_only=False, distance_threshold=50_000):
        self.add_data = with_data
        self.bool_only = bool_only
        self.data = list()
        self.last_r_start = -2_147_483_648
        self.last_q_start = -2_147_483_648
        self.current_r_end = 0
        self.index = 0
        self.distance_threshold = distance_threshold

    def __len__(self):
        return self.starts.size()

    def set_distance_threshold(self, threshold):
        self.distance_threshold = threshold

    cpdef add(self, int start, int end, value=None):
        """Add an interval. A value can also be accepted if with_data=True
        Input intervals must be added in sorted order
        
        :param start: The start of the interval
        :type start: int
        :param end: The start of the interval
        :type end: int
        :param value: Any associated value
        :type value: object
        """
        if end < start:
            raise ValueError("End less than start")
        if start < self.last_r_start:
            raise ValueError(f'Interval {start}-{end} is not in sorted order, last interval seen was {self.last_r_start}')
        self.starts.push_back(start)
        self.current_r_end = max(self.current_r_end, end)
        self.ends.push_back(Interval(end, self.current_r_end))
        self.last_r_start = start
        if self.add_data:
            self.data.append(value)

    cpdef add_from_iter(self, iterable):
        """Add intervals from an iterable. A value can also be accepted if with_data=True
        Input intervals must be added in sorted order
        
        :param iterable: Input intervals (start, end) or (start, end, value)
        :type iterable: object
        """
        cdef int start, end
        if hasattr(iterable, "len"):
            self.starts.reserve(len(iterable))
            self.ends.reserve(len(iterable))
        for item in iterable:
            start = item[0]
            end = item[1]
            assert end >= start
            if start < self.last_r_start:
                raise ValueError(f'Interval {start}-{end} is not in sorted order, last interval seen was {self.last_r_start}')
            self.starts.push_back(start)
            self.current_r_end = max(self.current_r_end, end)
            self.ends.push_back(Interval(end, self.current_r_end))
            self.last_r_start = start
            if self.add_data:
                self.data.append(item[2])

    def items(self):
        """Returns an interator over the stored intervals
        Yields:
            (start, end)        # with_data=False
            (start, end, data)  # with_data=True
        """
        cdef size_t i
        for i in range(self.starts.size()):
            if self.add_data:
                yield self.starts[i], self.ends[i].end, self.data[i]
            else:
                yield self.starts[i], self.ends[i].end

    cdef void _line_scan(self, int pos):
        cdef size_t i
        if pos < self.last_q_start:
            i = self.index
            while i > 0 and pos <= self.starts[i]:
                i -= 1
            while True:
                if self.ends[i].covered_end >= pos:
                    if i > 0:
                        i -= 1
                    else:
                        break
                else:
                    if i < self.ends.size() - 1:
                        i += 1
                    break
            self.index = i
        else:
            i = self.index
            while i < self.ends.size() and pos > self.ends[i].end:
                i += 1
            self.index = i

    cdef void _binary_search(self, int pos):
        cdef vector[int].iterator lower
        if self.last_q_start < pos:
            lower = lower_bound(self.starts.begin() + self.index, self.starts.end(), pos)
        else:
            lower = lower_bound(self.starts.begin(), self.starts.begin() + self.index, pos)
        self.index = lower - self.starts.begin()
        while True:
            if self.ends[self.index].covered_end >= pos:
                if self.index > 0:
                    self.index -= 1
                else:
                    break
            else:
                if self.index < self.ends.size() - 1:
                    self.index += 1
                break
    cdef void _set_reference_index(self, int pos):
        # self._binary_search(pos)
        # if self.query_sorted:
        #     self._line_scan(pos)
        # else:
        if c_abs(pos - self.last_q_start) > self.distance_threshold:
            self._binary_search(pos)
        else:
            self._line_scan(pos)
    cpdef search_point(self, int pos):
        """Search for reference intervals overlapping single query point

        :param pos: The position to search for overlaps
        :type pos: int
        """
        cdef size_t size = self.ends.size()
        cdef size_t i
        if size == 0:
            return False if self.bool_only else []
        self._set_reference_index(pos)
        self.last_q_start = pos
        if not self.bool_only:
            found = []
            for i in range(self.index, size):
                if self.starts[i] <= pos <= self.ends[i].end:
                    if self.add_data:
                        found.append((self.starts[i], self.ends[i].end, self.data[i]))
                    else:
                        found.append((self.starts[i], self.ends[i].end))
                    continue
                elif pos < self.starts[i] or pos > self.ends[i].covered_end:
                    break
            return found
        else:
            return False if self.index >= size else self.starts[self.index] <= pos <= self.ends[self.index].end

    cpdef search_interval(self, int start, int end):
        """Search for reference intervals overlapping single query interval

        :param start: The start position to search for overlaps
        :type start: int
        :param end: The end position to search for overlaps
        :type end: int
        """
        cdef size_t size = self.ends.size()
        cdef size_t i
        if size == 0:
            return False if self.bool_only else []
        self._set_reference_index(start)
        self.last_q_start = start
        if not self.bool_only:
            found = []
            for i in range(self.index, size):
                if is_overlapping(start, end, self.starts[i], self.ends[i].end):
                    if self.add_data:
                        found.append((self.starts[i], self.ends[i].end, self.data[i]))
                    else:
                        found.append((self.starts[i], self.ends[i].end))
                    continue
                elif start < self.starts[i] or start > self.ends[i].covered_end:
                    break
            return found
        else:
            return False if self.index >= size else is_overlapping(start, end, self.starts[self.index], self.ends[self.index].end)
