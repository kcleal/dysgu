// Version 0.3.7
#pragma once

#include <algorithm>
#include <vector>
#include <climits>
#include <cstdint>
#include <iostream>
#include <limits>
#ifndef SI_NOSIMD
    #if defined(__AVX2__)
        #include <immintrin.h>
    #elif defined(__ARM_NEON__) || defined(__aarch64__)
        #include <arm_neon.h>
    #else
        #define SI_NOSIMD
    #endif
#endif

namespace si {

// Core data structure
template<typename S, typename T>
struct Interval {
    S start, end;
    T data;
    Interval() = default;
    Interval(S s, S e, T d) : start(s), end(e), data(d) {}
};

/**
 * @file SuperIntervals.hpp
 * @brief IntervalMap is an associative data structure for finding interval intersections
 *
 * IntervalMap is a template class that provides efficient interval intersection operations.
 * It supports adding intervals, indexing them for fast queries, and performing various
 * intersection operations using an implicit superinterval tree algorithm.
 *
 * @note Intervals are considered end-inclusive
 * @note The build() function must be called before any queries. If more intervals are added, call build() again.
 *
 * @tparam S The scalar type for interval start and end points (e.g., int, float)
 * @tparam T The data type associated with each interval
 */
template<typename S, typename T>

class IntervalMap {
    public:
    std::vector<S> starts;
//    std::vector<S> ends;

    #ifdef __AVX2__
        alignas(32) std::vector<S> ends;
    #elif defined(__ARM_NEON)
        alignas(16) std::vector<S> ends;
    #else
        std::vector<S> ends;
    #endif
    std::vector<size_t> branch;
    std::vector<T> data;
    bool start_sorted, end_sorted;

    IntervalMap() : start_sorted(true), end_sorted(true) {}
    virtual ~IntervalMap() = default;

    /**
     * @brief Clears all intervals and resets the data structure
     */
    void clear() noexcept {
        data.clear(); starts.clear(); ends.clear(); branch.clear(); // idx = 0;
    }

    /**
     * @brief Reserves memory for a specified number of intervals
     * @param n Number of intervals to reserve space for
     */
    void reserve(size_t n) {
        data.reserve(n); starts.reserve(n); ends.reserve(n);
    }

    /**
    * @brief Returns the number of intervals in the data structure
    * @return Number of intervals
    */
    size_t size() {
        return starts.size();
    }

    /**
     * @brief Adds a new interval to the data structure
     * @param start Start point of the interval
     * @param end End point of the interval
     * @param value Data associated with the interval
     */
    void add(S start, S end, const T& value) {
        if (start_sorted && !starts.empty()) {
            start_sorted = (start < starts.back()) ? false : true;
            if (start_sorted && start == starts.back() && end > ends.back()) {
                end_sorted = false;
            }
        }
        starts.push_back(start);
        ends.push_back(end);
        data.emplace_back(value);
    }

    /**
     * @brief Build the index.
     * This function must be called after adding intervals and before performing any queries.
     * If more intervals are added after indexing, this function should be called again.
     */
    virtual void build() {
        if (starts.size() == 0) {
            return;
        }
        sort_intervals();
        branch.resize(starts.size(), SIZE_MAX);
        std::vector<std::pair<S, size_t>> br;
        br.reserve((starts.size() / 10) + 1);
        br.emplace_back() = {ends[0], 0};
        for (size_t i=1; i < ends.size(); ++i) {
            while (!br.empty() && br.back().first < ends[i]) {
                br.pop_back();
            }
            if (!br.empty()) {
                branch[i] = br.back().second;
            }
            br.emplace_back() = {ends[i], i};
        }
    }

    /**
     * @brief Retrieves an interval at a specific index
     * @param index The index of the interval to retrieve
     * @return The Interval at the specified index
     */
    const Interval<S, T>& at(size_t index) const {
        return Interval<S, T>{starts[index], ends[index], data[index]};
    }
    /**
     * @brief Retrieves an interval at a specific index and writes it into the provided Interval object.
     * @param index The index of the interval to retrieve.
     * @param itv  Reference to an Interval object to populate with the retrieved interval.
     */
    void at(size_t index, Interval<S, T>& itv) {
        itv.start = starts[index];
        itv.end = ends[index];
        itv.data = data[index];
    }

    // Iterator interfaces
    class IndexIterator {
    public:
        IndexIterator(const IntervalMap* parent, size_t pos, S query_start)
            : parent_(parent), current_pos_(pos), query_start_(query_start),
              has_value_(false), value_(0) {
            next();
        }

        size_t operator*() const { return value_; }

        IndexIterator& operator++() {
            next();
            return *this;
        }

        bool operator!=(const IndexIterator& other) const {
            return has_value_ != other.has_value_;
        }

        bool operator==(const IndexIterator& other) const {
            return has_value_ == other.has_value_;
        }

    private:
        void next() const {
            if (current_pos_ == SIZE_MAX) {
                has_value_ = false;
                return;
            }
            if (query_start_ <= parent_->ends[current_pos_]) {
                value_ = parent_->data[current_pos_];
                current_pos_ -= 1;
                has_value_ = true;
                return;
            }
            while (true) {
                current_pos_ =parent_->branch[current_pos_];
                if (current_pos_ == SIZE_MAX) {
                    break;
                }
                if (query_start_ <= parent_->ends[current_pos_]) {
                    value_ = parent_->data[current_pos_];
                    current_pos_ -= 1;
                    has_value_ = true;
                    return;
                }
            }
            has_value_ = false;
            return;
        }

        const IntervalMap* parent_;
        mutable size_t current_pos_;
        S query_start_;
        mutable bool has_value_;
        mutable size_t value_;
    };

    class IndexRange {
    public:
        IndexRange(const IntervalMap* parent, S start, S end)
            : parent_(parent), start_(start), end_(end) {}
        IndexIterator begin() const {
            const size_t pos = parent_->starts.empty() ? SIZE_MAX : parent_->upper_bound(end_);
            return IndexIterator(parent_, pos, start_);
        }
        IndexIterator end() const {
            return IndexIterator(parent_, SIZE_MAX, start_);
        }
    private:
        const IntervalMap* parent_;
        S start_, end_;
    };

    /**
     * @brief Creates a range object for iterating over interval indices that intersect [start, end]
     * @param start Start point of the search range
     * @param end End point of the search range
     * @return IndexRange object that can be used with range-based for loops
     */
    IndexRange search_idxs(S start, S end) const noexcept {
        return IndexRange(this, start, end);
    }

    class KeyIterator {
    public:
        KeyIterator(const IntervalMap* parent, size_t pos, S query_start)
            : parent_(parent), current_pos_(pos), query_start_(query_start),
              has_value_(false), value_{} {
            next();
        }

        const std::pair<S, S>& operator*() const { return value_; }

        KeyIterator& operator++() {
            next();
            return *this;
        }

        bool operator!=(const KeyIterator& other) const {
            return has_value_ != other.has_value_;
        }

        bool operator==(const KeyIterator& other) const {
            return has_value_ == other.has_value_;
        }

    private:
        void next() const {
            if (current_pos_ == SIZE_MAX) {
                has_value_ = false;
                return;
            }
            if (query_start_ <= parent_->ends[current_pos_]) {
                value_.first = parent_->starts[current_pos_];
                value_.second = parent_->ends[current_pos_];
                current_pos_ -= 1;
                has_value_ = true;
                return;
            }
            while (true) {
                current_pos_ =parent_->branch[current_pos_];
                if (current_pos_ == SIZE_MAX) {
                    break;
                }
                if (query_start_ <= parent_->ends[current_pos_]) {
                    value_.first = parent_->starts[current_pos_];
                    value_.second = parent_->ends[current_pos_];
                    current_pos_ -= 1;
                    has_value_ = true;
                    return;
                }
            }
            has_value_ = false;
            return;
        }

        const IntervalMap* parent_;
        mutable size_t current_pos_;
        S query_start_;
        mutable bool has_value_;
        mutable std::pair<S, S> value_;
    };

    class KeyRange {
    public:
        KeyRange(const IntervalMap* parent, S start, S end)
            : parent_(parent), start_(start), end_(end) {}
        KeyIterator begin() const {
            const size_t pos = parent_->starts.empty() ? SIZE_MAX : parent_->upper_bound(end_);
            return KeyIterator(parent_, pos, start_);
        }
        KeyIterator end() const {
            return KeyIterator(parent_, SIZE_MAX, start_);
        }
    private:
        const IntervalMap* parent_;
        S start_, end_;
    };

    /**
     * @brief Creates a range object for iterating over interval keys that intersect [start, end]
     * @param start Start point of the search range
     * @param end End point of the search range
     * @return KeyRange object that can be used with range-based for loops
     */
    KeyRange search_keys(S start, S end) const noexcept {
        return KeyRange(this, start, end);
    }


    class ValueIterator {
    public:
        ValueIterator(const IntervalMap* parent, size_t pos, S query_start)
            : parent_(parent), current_pos_(pos), query_start_(query_start),
              has_value_(false), value_{} {
            next();
        }

        const T& operator*() const { return value_; }

        ValueIterator& operator++() {
            next();
            return *this;
        }

        bool operator!=(const ValueIterator& other) const {
            return has_value_ != other.has_value_;
        }

        bool operator==(const ValueIterator& other) const {
            return has_value_ == other.has_value_;
        }

    private:
        void next() const {
            if (current_pos_ == SIZE_MAX) {
                has_value_ = false;
                return;
            }
            if (query_start_ <= parent_->ends[current_pos_]) {
                value_ = parent_->data[current_pos_];
                current_pos_ -= 1;
                has_value_ = true;
                return;
            }
            while (true) {
                current_pos_ =parent_->branch[current_pos_];
                if (current_pos_ == SIZE_MAX) {
                    break;
                }
                if (query_start_ <= parent_->ends[current_pos_]) {
                    value_ = parent_->data[current_pos_];
                    current_pos_ -= 1;
                    has_value_ = true;
                    return;
                }
            }
            has_value_ = false;
            return;
        }

        const IntervalMap* parent_;
        mutable size_t current_pos_;
        S query_start_;
        mutable bool has_value_;
        mutable T value_;
    };

    class ValueRange {
    public:
        ValueRange(const IntervalMap* parent, S start, S end)
            : parent_(parent), start_(start), end_(end) {}
        ValueIterator begin() const {
            const size_t pos = parent_->starts.empty() ? SIZE_MAX : parent_->upper_bound(end_);
            return ValueIterator(parent_, pos, start_);
        }
        ValueIterator end() const {
            return ValueIterator(parent_, SIZE_MAX, start_);
        }
    private:
        const IntervalMap* parent_;
        S start_, end_;
    };

    /**
     * @brief Creates a range object for iterating over interval values that intersect [start, end]
     * @param start Start point of the search range
     * @param end End point of the search range
     * @return ValueRange object that can be used with range-based for loops
     */
    ValueRange search_values(S start, S end) const noexcept {
        return ValueRange(this, start, end);
    }

    class ItemIterator {
    public:
        ItemIterator(const IntervalMap* parent, size_t pos, S query_start)
            : parent_(parent), current_pos_(pos), query_start_(query_start),
              has_value_(false), value_{} {
            next();
        }

        const Interval<S, T>& operator*() const { return value_; }

        ItemIterator& operator++() {
            next();
            return *this;
        }

        bool operator!=(const ItemIterator& other) const {
            return has_value_ != other.has_value_;
        }

        bool operator==(const ItemIterator& other) const {
            return has_value_ == other.has_value_;
        }

    private:
        void next() const {
            if (current_pos_ == SIZE_MAX) {
                has_value_ = false;
                return;
            }
            if (query_start_ <= parent_->ends[current_pos_]) {
                value_.start = parent_->starts[current_pos_];
                value_.end = parent_->ends[current_pos_];
                value_.data = parent_->data[current_pos_];
                current_pos_ -= 1;
                has_value_ = true;
                return;
            }
            while (true) {
                current_pos_ =parent_->branch[current_pos_];
                if (current_pos_ == SIZE_MAX) {
                    break;
                }
                if (query_start_ <= parent_->ends[current_pos_]) {
                    value_.start = parent_->starts[current_pos_];
                    value_.end = parent_->ends[current_pos_];
                    value_.data = parent_->data[current_pos_];
                    current_pos_ -= 1;
                    has_value_ = true;
                    return;
                }
            }
            has_value_ = false;
            return;
        }

        const IntervalMap* parent_;
        mutable size_t current_pos_;
        S query_start_;
        mutable bool has_value_;
        mutable Interval<S, T> value_;
    };

    class ItemRange {
    public:
        ItemRange(const IntervalMap* parent, S start, S end)
            : parent_(parent), start_(start), end_(end) {}
        ItemIterator begin() const {
            const size_t pos = parent_->starts.empty() ? SIZE_MAX : parent_->upper_bound(end_);
            return ItemIterator(parent_, pos, start_);
        }
        ItemIterator end() const {
            return ItemIterator(parent_, SIZE_MAX, start_);
        }
    private:
        const IntervalMap* parent_;
        S start_, end_;
    };

    /**
     * @brief Creates a range object for iterating over interval items that intersect [start, end]
     * @param start Start point of the search range
     * @param end End point of the search range
     * @return ItemRange object that can be used with range-based for loops
     */
    ItemRange search_items(S start, S end) const noexcept {
        return ItemRange(this, start, end);
    }

    /**
     * @brief Finds the largest index such that starts[index] ≤ value.
     * @param value The upper bound value to search for.
     * @note After calling this, idx will be set to that index (or SIZE_MAX if none).
     */
    virtual inline size_t upper_bound(const S value) const noexcept {
        size_t length = starts.size();
        size_t idx = 0;
        while (length > 1) {
            const size_t half = length / 2;
            idx += (starts[idx + half] <= value) * (length - half);
            length = half;
        }
        if (starts[idx] > value) {
            --idx;  // Might underflow to SIZE_MAX
        }
        return idx;

//        idx = std::distance(starts.begin(),
//            std::upper_bound(starts.begin(), starts.end(), value)) - 1;
    }

    /**
     * @brief Narrows a [left…right) range so that starts[left] is the first element ≥ value.
     * @param value The lower‐bound value to search for.
     * @param left  On entry, the lower end of the search range; on exit, the found index or SIZE_MAX.
     * @param right The upper end of the search range (exclusive).
     */
    inline void upper_bound_range(const S value, size_t& left, const size_t right) const noexcept {
        // First do exponential search if we have room
        size_t search_right = right;
        size_t bound = 1;
        while (left > 0 && value < starts[left]) {
            search_right = left;
            left = (bound <= left) ? left - bound : 0;
            bound *= 2;
        }
        // Now do binary search in the range
        size_t length = search_right - left;
        while (length > 1) {
            const size_t half = length / 2;
            left += (starts[left + half] < value) * (length - half);
            length = half;
        }
        if (left == 0 && starts[left] >= value) {
            left = SIZE_MAX;
        }
    }

    /**
     * @brief Collects all data values whose intervals intersect [start,end].
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @param found Output vector that will be filled with matching data values.
     */
    void search_values(const S start, const S end, std::vector<T>& found) const {
        if (starts.empty()) {
            return;
        }
        const size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return;
        }
        size_t i = idx;
        while (i != SIZE_MAX && start <= ends[i]) {
            --i;
        }
        if (i == SIZE_MAX) {
            found.insert(found.end(), data.rend() - idx - 1, data.rend());
            return;
        }
        found.insert(found.end(), data.rend() - idx - 1, data.rend() - i - 1);
        // This is slightly faster, but would lose the sort order:
        // found.insert(found.end(), data.begin() + i, data.begin() + idx);
        i = branch[i];
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                found.push_back(data[i]);
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    /**
     * @brief Like search_intervals, but optimized for large query ranges.
     * @note Uses exponential search - best when query range is large relative to stored intervals
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @param found Output vector that will be filled with matching data values.
     */
    void search_values_large(const S start, const S end, std::vector<T>& found) const {
        if (starts.empty()) {
            return;
        }
        const size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return;
        }
        size_t i = idx;
        upper_bound_range(start, i, idx);
        while (i != SIZE_MAX && start <= ends[i]) {
            --i;
        }
        if (i == SIZE_MAX) {
            found.insert(found.end(), data.rend() - idx - 1, data.rend());
            return;
        }
        found.insert(found.end(), data.rend() - idx - 1, data.rend() - i - 1);
        i = branch[i];
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                found.push_back(data[i]);
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    /**
     * @brief Counts how many intervals intersect [start,end].
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @return Number of intervals overlapping the query.
     */
    size_t count_linear(const S start, const S end) const noexcept {
        if (starts.empty()) {
            return 0;
        }
        size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return 0;
        }
        size_t starting_i = idx;
        while (idx != SIZE_MAX && start <= ends[idx]) {
            --idx;
        }
        size_t count = starting_i - idx;
        if (idx == SIZE_MAX) {
            return count;
        }
        idx = branch[idx];
        while (idx != SIZE_MAX) {
            if (start <= ends[idx]) {
                count += 1;
                --idx;
            } else {
                idx = branch[idx];
            }
        }
        return count;
    }

    size_t count(const S start, const S end) const noexcept {
        if (starts.empty()) {
            return 0;
        }
        size_t i = upper_bound(end);
        if (i == SIZE_MAX) {
            return 0;
        }
        size_t found = 0;

        #if defined(SI_NOSIMD)

        constexpr size_t block = 16;
        constexpr simd_kind active_simd = simd_kind::none;

        #elif defined(__AVX2__)

        __m256i start_vec = _mm256_set1_epi32(start);
        constexpr size_t simd_width = 256 / (sizeof(S) * 8);
        constexpr size_t block = simd_width * 4;  // 2 cache lines
        constexpr simd_kind active_simd = simd_kind::avx2;

        #elif defined(__ARM_NEON__) || defined(__aarch64__)

        int32x4_t start_vec = vdupq_n_s32(start);
        constexpr size_t simd_width = 128 / (sizeof(S) * 8);
        uint32x4_t ones = vdupq_n_u32(1);
        constexpr size_t block = simd_width * 8;  // 2 cache lines
        constexpr simd_kind active_simd = simd_kind::neon;

        #endif

        while (i > 0) {
            if (start <= ends[i]) {
                ++found;
                --i;
                // Types with width !=4 will use the no-simd path here
                if constexpr (active_simd == simd_kind::none || sizeof(S) != 4) {
                    while (i > block) {
                        size_t count = 0;
                        for (size_t j = i; j > i - block; --j) {
                            count += (start <= ends[j]) ? 1 : 0;
                        }
                        found += count;
                        i -= block;
                        if (count < block && start > ends[i + 1]) {  // check for a branch
                            break;
                        }
                    }
                }
                #if defined(__AVX2__)
                else if constexpr (active_simd == simd_kind::avx2) {
//                    while (i > block) {
//                        size_t count = 0;
//                        for (size_t j = i - block + 1; j < i; j += simd_width) {
//                            __m256i ends_vec = _mm256_loadu_si256((__m256i*)(&ends[j - simd_width + 1]));
//                            __m256i cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
//                            int mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask));
//                            count += 8 - _mm_popcnt_u32(mask);
//                        }
//                        found += count;
//                        i -= block;
//                        if (count < block) {
//                            break;
//                        }
//                    }

                    while (i > block) {
                        size_t j = i - block + 1;

                        // Load all 4 vectors
                        __m256i ends_vec0 = _mm256_load_si256((__m256i*)(&ends[j]));
                        __m256i ends_vec1 = _mm256_load_si256((__m256i*)(&ends[j + simd_width]));
                        __m256i ends_vec2 = _mm256_load_si256((__m256i*)(&ends[j + 2 * simd_width]));
                        __m256i ends_vec3 = _mm256_load_si256((__m256i*)(&ends[j + 3 * simd_width]));

                        // Compare all vectors
                        __m256i cmp_mask0 = _mm256_cmpgt_epi32(start_vec, ends_vec0);
                        __m256i cmp_mask1 = _mm256_cmpgt_epi32(start_vec, ends_vec1);
                        __m256i cmp_mask2 = _mm256_cmpgt_epi32(start_vec, ends_vec2);
                        __m256i cmp_mask3 = _mm256_cmpgt_epi32(start_vec, ends_vec3);

                        // Extract masks
                        int mask0 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask0));
                        int mask1 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask1));
                        int mask2 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask2));
                        int mask3 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask3));

                        // Count and accumulate
                        size_t count = (8 - _mm_popcnt_u32(mask0)) + (8 - _mm_popcnt_u32(mask1)) +
                                       (8 - _mm_popcnt_u32(mask2)) + (8 - _mm_popcnt_u32(mask3));

                        found += count;
                        i -= block;
                        if (count < block) {
                            break;
                        }
                    }
                }
                #elif defined(__ARM_NEON__) || defined(__aarch64__)
                else {  // NEON
//                    while (i > block) {
//                        size_t count = 0;
//                        uint32x4_t mask, bool_mask;
//                        for (size_t j = i - block + 1; j < i; j += simd_width) { // Neon 4 int32 at a time
//                            int32x4_t ends_vec = vld1q_s32(&ends[j]);
//                            mask = vcgtq_s32(start_vec, ends_vec);  // start > ends[j]
//                            bool_mask = vaddq_u32(mask, ones);
//                            count += vaddvq_u32(bool_mask); // Sum all lanes
//                        }
//                        found += count;
//                        i -= block;
//                        if (count < block) {  // check for overlap again, before checking for branch?
//                            break;
//                        }
//                    }
                    while (i > block) {
                        size_t j = i - block + 1;

                        // Load all 8 vectors
                        int32x4_t ends_vec0 = vld1q_s32(&ends[j]);
                        int32x4_t ends_vec1 = vld1q_s32(&ends[j + simd_width]);
                        int32x4_t ends_vec2 = vld1q_s32(&ends[j + 2 * simd_width]);
                        int32x4_t ends_vec3 = vld1q_s32(&ends[j + 3 * simd_width]);
                        int32x4_t ends_vec4 = vld1q_s32(&ends[j + 4 * simd_width]);
                        int32x4_t ends_vec5 = vld1q_s32(&ends[j + 5 * simd_width]);
                        int32x4_t ends_vec6 = vld1q_s32(&ends[j + 6 * simd_width]);
                        int32x4_t ends_vec7 = vld1q_s32(&ends[j + 7 * simd_width]);

                        // Compare all vectors
                        uint32x4_t mask0 = vcgtq_s32(start_vec, ends_vec0);
                        uint32x4_t mask1 = vcgtq_s32(start_vec, ends_vec1);
                        uint32x4_t mask2 = vcgtq_s32(start_vec, ends_vec2);
                        uint32x4_t mask3 = vcgtq_s32(start_vec, ends_vec3);
                        uint32x4_t mask4 = vcgtq_s32(start_vec, ends_vec4);
                        uint32x4_t mask5 = vcgtq_s32(start_vec, ends_vec5);
                        uint32x4_t mask6 = vcgtq_s32(start_vec, ends_vec6);
                        uint32x4_t mask7 = vcgtq_s32(start_vec, ends_vec7);

                        // Convert to boolean masks
                        uint32x4_t bool_mask0 = vaddq_u32(mask0, ones);
                        uint32x4_t bool_mask1 = vaddq_u32(mask1, ones);
                        uint32x4_t bool_mask2 = vaddq_u32(mask2, ones);
                        uint32x4_t bool_mask3 = vaddq_u32(mask3, ones);
                        uint32x4_t bool_mask4 = vaddq_u32(mask4, ones);
                        uint32x4_t bool_mask5 = vaddq_u32(mask5, ones);
                        uint32x4_t bool_mask6 = vaddq_u32(mask6, ones);
                        uint32x4_t bool_mask7 = vaddq_u32(mask7, ones);

                        // Sum all lanes and accumulate
                        size_t count = vaddvq_u32(bool_mask0) + vaddvq_u32(bool_mask1) +
                                       vaddvq_u32(bool_mask2) + vaddvq_u32(bool_mask3) +
                                       vaddvq_u32(bool_mask4) + vaddvq_u32(bool_mask5) +
                                       vaddvq_u32(bool_mask6) + vaddvq_u32(bool_mask7);

                        found += count;
                        i -= block;
                        if (count < block) {
                            break;
                        }
                    }
                }
                #endif
            } else {
                if (branch[i] == SIZE_MAX) {
                    return found;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0] && starts[0] <= end) {
            ++found;
        }
        return found;
    }

    /**
     * @brief Like count, but optimized for large query ranges.
     * @note Uses exponential search - best when query range is large relative to stored intervals
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @return Number of intervals overlapping the query.
     */
    size_t count_large(const S start, const S end) const noexcept {
        if (starts.empty()) {
            return 0;
        }
        const size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return 0;
        }
        size_t i = idx;
        upper_bound_range(start, i, idx);
        while (i != SIZE_MAX && start <= ends[i]) {
            --i;
        }
        size_t count = idx - i;
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                count += 1;
                --i;
            } else {
                i = branch[i];
            }
        }
        return count;
    }

    /**
     * @brief Tests whether any interval overlaps the point range [start,end].
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @return true if at least one overlap exists, false otherwise.
     */
    bool has_overlaps(const S start, const S end) const noexcept {
        if (starts.empty()) {
            return false;
        }
        const size_t idx = upper_bound(end);
        return idx != SIZE_MAX && start <= ends[idx];
    }

    /**
     * @brief Collects the indices of all intervals intersecting [start,end].
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @param found Output vector that will be filled with matching interval indices.
     */
    void search_idxs(const S start, const S end, std::vector<size_t>& found) const {
        if (starts.empty()) {
            return;
        }
        const size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return;
        }
        size_t i = idx;
        while (i != SIZE_MAX && start <= ends[i]) {
            --i;
        }
        if (i == SIZE_MAX) {
            found.insert(found.end(), CountingIterator(0), CountingIterator(idx + 1));
            return;
        }
        found.insert(found.end(), CountingIterator(i + 1), CountingIterator(idx + 1));
        i = branch[i];
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                found.push_back(i);
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    /**
     * @brief Collects the (start,end) pairs of all intervals intersecting [start,end].
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @param found Output vector that will be filled with matching interval key pairs.
     */
    void search_keys(const S start, const S end, std::vector<std::pair<S, S>>& found) const {
        if (starts.empty()) {
            return;
        }
        const size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return;
        }
        size_t i = idx;
        while (i != SIZE_MAX && start <= ends[i]) {
            found.emplace_back() = {starts[i], ends[i]};
            --i;
        }
        if (i == SIZE_MAX) {
            return;
        }
        i = branch[i];
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                found.emplace_back() = {starts[i], ends[i]};
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    /**
     * @brief Collects the full Interval objects intersecting [start,end].
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @param found Output vector that will be filled with matching Interval instances.
     */
    void search_items(const S start, const S end, std::vector<Interval<S, T>>& found) const {
        if (starts.empty()) {
            return;
        }
        const size_t idx = upper_bound(end);
        if (idx == SIZE_MAX) {
            return;
        }
        size_t i = idx;
        while (i != SIZE_MAX && start <= ends[i]) {
            found.emplace_back() = {starts[i], ends[i], data[i]};
            --i;
        }
        if (i == SIZE_MAX) {
            return;
        }
        i = branch[i];
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                found.emplace_back() = {starts[i], ends[i], data[i]};
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    /**
     * @brief Computes how many intervals overlap [start,end] and the total covered length.
     * @param start      The start of the query interval.
     * @param end        The end of the query interval.
     * @param cov_result Pair where first = count of overlaps, second = sum of overlapping lengths.
     */
    void coverage(const S start, const S end, std::pair<size_t, S> &cov_result) const {
        if (starts.empty()) {
            return;
        }
        size_t i = upper_bound(end);
        if (i == SIZE_MAX) {
            return;
        }
//        size_t i = idx;
        while (i != SIZE_MAX && start <= ends[i]) {
            ++cov_result.first;
            cov_result.second += std::min(ends[i], end) - std::max(starts[i], start);
            --i;
        }
        if (i == SIZE_MAX) {
            return;
        }
        i = branch[i];
        while (i != SIZE_MAX) {
            if (start <= ends[i]) {
                ++cov_result.first;
                cov_result.second += std::min(ends[i], end) - std::max(starts[i], start);
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    /**
     * @brief Finds all data values of intervals “stabbed” by a single point.
     * @param point The point to test.
     * @param found Output vector that will be filled with data values of intervals containing point.
     */
    void search_point(const S point, std::vector<T>& found) const {
        if (starts.empty()) {
            return;
        }
        size_t i = upper_bound(point);
//        size_t i = idx;
        while (i != SIZE_MAX && point <= ends[i]) {
            found.push_back(data[i]);
            --i;
        }
        i = branch[i];
        while (i != SIZE_MAX) {
            if (point <= ends[i]) {
                found.push_back(data[i]);
                --i;
            } else {
                i = branch[i];
            }
        }
    }

    protected:

    enum class simd_kind { none, avx2, neon };

    std::vector<Interval<S, T>> tmp;

    template<typename CompareFunc>
    void sort_block(size_t start_i, size_t end_i, CompareFunc compare) {
        size_t range_size = end_i - start_i;
        tmp.resize(range_size);
        for (size_t i = 0; i < range_size; ++i) {
            tmp[i] = Interval<S, T>(starts[start_i + i], ends[start_i + i], data[start_i + i]);
        }
        std::sort(tmp.begin(), tmp.end(), compare);
        for (size_t i = 0; i < range_size; ++i) {
            starts[start_i + i] = tmp[i].start;
            ends[start_i + i] = tmp[i].end;
            data[start_i + i] = tmp[i].data;
        }
    }

    /**
     * @brief Ensures the global interval list is properly sorted by start (and end within ties).
     *        If only ends need resorting among equal starts, does a more targeted pass.
     */
    void sort_intervals() {
        if (!start_sorted) {
            sort_block(0, starts.size(),
                [](const Interval<S, T>& a, const Interval<S, T>& b) { return (a.start < b.start || (a.start == b.start && a.end > b.end)); });
            start_sorted = true;
            end_sorted = true;
        } else if (!end_sorted) {  // only sort parts that need sorting - ends in descending order
            size_t it_start = 0;
            while (it_start < starts.size()) {
                size_t block_end = it_start + 1;
                bool needs_sort = false;
                while (block_end < starts.size() && starts[block_end] == starts[it_start]) {
                    if (block_end > it_start && ends[block_end] > ends[block_end - 1]) {
                        needs_sort = true;
                    }
                    ++block_end;
                }
                if (needs_sort) {
                    sort_block(it_start, block_end, [](const Interval<S, T>& a, const Interval<S, T>& b) { return a.end > b.end; });
                }
                it_start = block_end;
            }
            end_sorted = true;
        }
    }

    struct CountingIterator {  // Just for using insert indexes
        size_t value;
        using iterator_category = std::forward_iterator_tag;
        using value_type = size_t;
        using difference_type = std::ptrdiff_t;
        using pointer = const size_t*;
        using reference = size_t;

        explicit CountingIterator(size_t v) : value(v) {}
        size_t operator*() const { return value; }
        CountingIterator& operator++() { ++value; return *this; }
        bool operator!=(const CountingIterator& other) const { return value != other.value; }
    };
};


template<typename S, typename T>
class IntervalMapEytz : public IntervalMap<S, T> {
public:
    /**
     * @brief Builds the Eytzinger‐layout index and branch pointers for fast searching.
     * @note Overrides the base build() to use the cache‐friendly layout.
     */
    void build() override {
        if (this->starts.size() == 0) {
            return;
        }
        this->starts.shrink_to_fit();
        this->ends.shrink_to_fit();
        this->data.shrink_to_fit();
        this->sort_intervals();

        eytz.resize(this->starts.size() + 1);
        eytz_index.resize(this->starts.size() + 1);
        eytzinger(&this->starts[0], this->starts.size());

        this->branch.resize(this->starts.size(), SIZE_MAX);
        std::vector<std::pair<S, size_t>> br;
        br.reserve(1000);
        br.emplace_back() = {this->ends[0], 0};
        for (size_t i=1; i < this->ends.size(); ++i) {
            while (!br.empty() && br.back().first < this->ends[i]) {
                br.pop_back();
            }
            if (!br.empty()) {
                this->branch[i] = br.back().second;
            }
            br.emplace_back() = {this->ends[i], i};
        }
//        this->idx = 0;
    }

    /**
     * @brief Uses the Eytzinger‐layout tree to locate the largest start ≤ x.
     * @param x The search value.
     * @note Overrides the base upper_bound to navigate the implicit binary tree.
     */
    inline size_t upper_bound(const S x) const noexcept override {
         size_t i = 0;
         const size_t n_intervals = this->starts.size();
         while (i < n_intervals) {
             if (eytz[i] > x) {
                 i = 2 * i + 1;
             } else {
                 i = 2 * i + 2;
             }
         }
         int shift = __builtin_ffs(~(i + 1));
         size_t best_idx = (i >> shift) - ((shift > 1) ? 1 : 0);
         i = (best_idx < n_intervals) ? eytz_index[best_idx] : n_intervals - 1;
         if (i > 0 && this->starts[i] > x) {
             --i;
         }
         return i;
    }

private:
    std::vector<S> eytz;
    std::vector<size_t> eytz_index;

    size_t eytzinger_helper(S* arr, size_t n, size_t i, size_t k) {
        if (k < n) {
            i = eytzinger_helper(arr, n, i, 2*k+1);
            eytz[k] = this->starts[i];
            eytz_index[k] = i;
            ++i;
            i = eytzinger_helper(arr, n, i, 2*k + 2);
        }
        return i;
    }

    int eytzinger(S* arr, size_t n) {
        return eytzinger_helper(arr, n, 0, 0);
    }
};

} // si
