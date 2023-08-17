#include <cstdint>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
#include <queue>
#include <map>
#include <cassert>
#include <cmath>
#include <vector>

#include "robin_hood.h"
#include "IITreeBFS.h"

typedef std::pair<int, uint8_t> PairW;  // v, weight
typedef std::pair<int, int> PairW2;  // v, weight


class DiGraph {
    // Directed weighted
    public:

        DiGraph() {}  // Construct
        ~DiGraph() {}

        std::vector<std::vector<PairW2>> inList;
        std::vector<std::vector<PairW2>> outList;

        int N = 0;
        int n_edges = 0;

        int addNode() {
            int n = N;
            std::vector<PairW2> node;
            std::vector<PairW2> in_node;
            outList.push_back(node);  // empty, no edges
            inList.push_back(in_node);
            N++;
            return n;
        }

        int hasEdge(int u, int v) {
            if ((u > N ) || (v > N)) { return 0; }
            for (const auto& val: outList[u]) {
                if (val.first == v) {
                    return 1;
                }
            }
            return 0;
        }

        void addEdge(int u, int v, int w) {
            // Assume hasEdge has been call
            outList[u].push_back(std::make_pair(v, w));
            inList[v].push_back(std::make_pair(u, w));
            n_edges += 1;
        }

        void updateEdge(int u, int v, int w) {
            // Check for edge and update, else create new edge
            if ((u > N ) || (v > N)) {
                outList[u].push_back(std::make_pair(v, w));
                inList[v].push_back(std::make_pair(u, w));
                n_edges += 1;
                return;
            }
            bool updated = false;
            for (auto& val: outList[u]) {
                if (val.first == v) {
                    val.second += w;
                    updated = true;
                }
            }
            if (updated == false) {
                outList[u].push_back(std::make_pair(v, w));
                n_edges += 1;
            }
            updated = false;
            for (auto& vali: inList[v]) {
                if (vali.first == u) {
                    vali.second += w;
                    updated = true;
                }
            }
            if (updated == false) {
                inList[v].push_back(std::make_pair(u, w));
            }
        }

        int weight(int u, int v) {
            if ((u > N ) || (v > N)) { return 0; }
            for (const auto& val: outList[u]) {
                if (val.first == v) {
                    return val.second;
                }
            }
            return 0;
        }

        int numberOfNodes() { return N; }

        void neighbors(int u, std::vector<int> &neigh) {  // Get all out edges
            //std::vector<int> neigh;
            if (!neigh.empty()) {
                neigh.clear();
            }
            for (const auto& val: outList[u]) {
                neigh.push_back(val.first);
            }
            //return neigh;
        }


        void forInEdgesOf(int u, std::vector<PairW2>& inEdges) {
//            std::vector<PairW2> inEdges;
            if (!inEdges.empty()) {
                inEdges.clear();
            }
            for (const auto& val: inList[u]) {
                inEdges.push_back(val);
            }
//            return inEdges;
        }

        float node_path_quality(int u, int v, int w) {
            int sum_in = 0;
            float q_in = 1;
            if (u != -1) {
                float in_weight = (float)weight(u, v);
                for (const auto& val: inList[v]) {
                    sum_in += val.second;
                }
                float other_in = sum_in - in_weight;
                if (other_in > 0) {
                    q_in = in_weight / sum_in;
                }
            }
            int sum_out = 0;
            float q_out = 1;
            if (w != -1) {
                float out_weight = (float)weight(v, w);
                for (const auto& val: outList[v]) {
                    sum_out += val.second;
                }
                float other_out = sum_out - out_weight;
                if (other_out > 0) {
                    q_out = out_weight / sum_out;
                }
            }
            if (q_out < q_in) { return q_out; };
            return q_in;
        }
};


class SimpleGraph {
    // Undirected, weighted
    public:

        SimpleGraph() {}  // Construct
        ~SimpleGraph() {}

        std::vector<std::vector<PairW>> adjList;

        int N = 0;
        int n_edges = 0;

        int addNode() {
            int n = N;
            std::vector<PairW> node;
            adjList.push_back(node);
            N++;
            return n;
        }

        int hasEdge(int u, int v) {
            if ((u > N ) || (v > N)) { return 0; }
            for (const auto& val: adjList[u]) {
                if (val.first == v) {
                    return 1;
                }
            }
            return 0;
        }

        void addEdge(int u, int v, uint8_t w) {
            adjList[u].push_back(std::make_pair(v, w));
            adjList[v].push_back(std::make_pair(u, w));
            n_edges += 1;
        }

        int edgeCount() { return n_edges; }

        int weight(int u, int v) {
            if ((u > N ) || (v > N)) { return 0; }
            for (const auto& val: adjList[u]) {
                if (val.first == v) {
                    return val.second;
                }
            }
            return 0;
        }

        void neighbors(int u, std::vector<int>& neigh) {
//            std::vector<int> neigh;
            if (!neigh.empty()) {
                neigh.clear();
            }
            for (const auto& val: adjList[u]) {
                neigh.push_back(val.first);
            }
//            return neigh;
        }

        void removeNode(int u) {
            // make dead nodes, takes more memory but faster
            for (auto& val_u: adjList[u]) {
                for (auto& val_v: adjList[val_u.first]) {
                    if (val_v.first == u) {
                        val_v.first = -1;
                        val_v.second = -1;
                    }
                }
            }
            // Set node edges to empty
            std::vector<PairW> node;
            node.push_back(std::make_pair(-1, -1));
            adjList[u] = node;
        }

        void connectedComponents(const char* outpath, bool low_mem, std::vector<int>& components) {
            std::string outpath_string = outpath;
            std::ofstream outf(outpath);
            std::vector<bool> visited(adjList.size(), false);
            if (!components.empty()) {
                components.clear();
            }
//            std::vector<int> components;
            for (int u=0; u<N; u++) {
                if (visited[u] == false) {
                    if (!low_mem) {
                        components.push_back(u);
                    } else {
                        outf.write((char*)&u, sizeof(int32_t));
                    }
                    visited[u] = true;
                    std::vector<int> queue;
                    // visit direct neighbors
                    if (adjList[u].size() > 0) {
                        for (const auto& val: adjList[u]) {
                            if ((val.first != -1) && (visited[val.first] == false)) {
                                queue.push_back(val.first);
                            }
                        }
                    }
                    while (!queue.empty()) {
                        int v = queue.back();
                        queue.pop_back();
                        if (visited[v] == false) {
                            if (!low_mem) {
                                components.push_back(v);
                            } else {
                                outf.write((char*)&v, sizeof(int32_t));
                            }
                            visited[v] = true;
                            for (const auto& val2: adjList[v]) {
                                if ((val2.first != -1) && (visited[val2.first] == false)) {
                                    queue.push_back(val2.first);
                                }
                            }
                        }
                    }
                    // -1 is end of component
                    if (!low_mem) {
                        components.push_back(-1);
                    } else {
                        int v = -1;
                        outf.write((char*)&v, sizeof(int32_t));
                    }
                }
            }
            outf.close();
//            return components;
        }

        std::size_t showSize() {
            size_t tot = sizeof(adjList);  // outer size
            for (const auto &v: adjList) {
                tot += sizeof(v);
                for (PairW v2: v) {
                    tot += sizeof(v2);
                }
            }
            return tot;  // bytes
        }
};


typedef std::pair<int, int> lookup_result;


class Int2IntMap
{
    public:

        Int2IntMap() {}
        ~Int2IntMap() {}

        int size() { return map.size(); }

        void insert(int key, int value) { map[key] = value; }
        void erase(int key) { map.erase(key); }

        int get(int key) { return map.at(key); }

        int has_key(int key) {
            if (map.find(key) == map.end()) { return 0; }
            else { return 1; }
        }

        lookup_result get_value(int key) {

            if (map.find(key) != map.end()) {
                return std::make_pair(1, map[key]);
            }
            else {
                return std::make_pair(0, 0);
            }
        }

    private:
        robin_hood::unordered_flat_map<int, int> map;
};


class IntSet
{
    public:

        IntSet() {}
        ~IntSet() {}

        int size() { return set.size(); }

        void insert(int key) { set.insert(key); }
        void erase(int key) { set.erase(key); }

        int has_key(int key) {
            if (set.find(key) == set.end()) { return 0; }
            else { return 1; }
        }

    private:
        robin_hood::unordered_set<int> set;
};


class TwoWayMap
{
    public:

        TwoWayMap() {}
        ~TwoWayMap() {}

        uint64_t key_2_64(char seq, uint64_t current_pos, uint64_t offset, uint64_t code) {
            uint64_t packed_data = 0;
            packed_data |= seq;  // Max 4 bits set per htslib
            packed_data |= current_pos << 4;  // 32 bits
            packed_data |= offset << 36;  // 23 bits
            packed_data |= code << 59;  // 4 bits
            return packed_data;
        }

        void insert_tuple_key(uint64_t packed_data, int index) {
            string_key[packed_data] = index;
            index_key_map.push_back(packed_data);
        }

        int has_tuple_key(uint64_t packed_data) {
            robin_hood::unordered_flat_map<uint64_t, int>::const_iterator got = string_key.find(packed_data);
            if (got == string_key.end()) {
                return 0;
            } else {
                last_key = got->first;
                last_index = got->second;
                return 1;
            }
        }

        int get_index_prev() { return last_index; };
        int get_key_prev() { return last_key; };

        void key_2_vec(uint64_t packed_data, std::vector<int>& v) {
            v[0] = packed_data & 15;  // 4 bits set to 1
            v[1] = (packed_data >> 4) & 4294967295;  // 1 x 32 bits
            v[2] = (packed_data >> 36) & 8388607;  // 1 x 23 bits
            v[3] = (packed_data >> 59) & 15;
        }

        void idx_2_vec(int index, std::vector<int>& v) {
            uint64_t packed_data = index_key_map[index];
            v[0] = packed_data & 15;  // 4 bits set to 1
            v[1] = (packed_data >> 4) & 4294967295;  // 1 x 32 bits
            v[2] = (packed_data >> 36) & 8388607;  // 1 x 23 bits
            v[3] = (packed_data >> 59) & 15;
        }

    private:
        uint64_t last_key = 0;  // Set this on call to has_tuple_key to prevent calculating twice
        int last_index = 0;
        robin_hood::unordered_flat_map<uint64_t, int> string_key;
        std::vector<uint64_t> index_key_map;

};


void graph_node_2_vec(uint64_t packed_data, std::vector<int>& v) {  // used for debugging
    // base, current_pos, offset, soft-clip-side
    v.push_back(packed_data & 15);
    v.push_back((packed_data >> 4) & 4294967295);
    v.push_back((packed_data >> 36) & 8388607);
    v.push_back((packed_data >> 59) & 15);
}


typedef robin_hood::unordered_set<long> set_of_long_t;
typedef robin_hood::unordered_map<long, set_of_long_t> mm_map_t;

class MinimizerTable
// Key's are read-names (integer), the table value is a set with each item a minimizer (long)
{
    public:
        MinimizerTable() {}
        ~MinimizerTable() {}

        int size() { return mm_map.size(); }

        void insert(long key, long value) { mm_map[key].insert(value); }

        void erase(long key) { mm_map.erase(key); }

        void erase_lower(long key, long value) { mm_map[key].erase(value); }

        int has_key(long key) {
            itr = mm_map.find(key);
            if (itr == mm_map.end()) { return 0; }
            else { return 1; }
        }

        set_of_long_t::iterator get_iterator_begin () {
            return itr->second.begin();
        }

        set_of_long_t::iterator get_iterator_end () {
            return itr->second.end();
        }

        int has_lower_key(long key2) {
            set_itr = itr->second.find(key2);
            if (set_itr == itr->second.end()) { return 0; }
            else { return 1; }
        }

        long get_lower() {
            return *set_itr;
        }

        set_of_long_t::iterator get_iterator() {
            return set_itr;
        }

    private:
        mm_map_t mm_map;
        mm_map_t::iterator itr;
        set_of_long_t::iterator set_itr;
};


struct Interval { int low, high; };

class BasicIntervalTree
// https://www.geeksforgeeks.org/interval-tree/
{
    public:

        BasicIntervalTree() {}
        ~BasicIntervalTree() {}

        void add(int start, int end, int index) {
            tree.add(start, end, index);
        }

        void index() {
            tree.index();
        }

        bool searchInterval(int pos, int pos2) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            if (a.size() == 0)
                return false;
            else
                return true;
        }

        void allOverlappingIntervals(int pos, int pos2, std::vector<int>& res) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            for (size_t i = 0; i < a.size(); ++i) {
                res.push_back(tree.data(a[i]));
            };
        }

        int countOverlappingIntervals(int pos, int pos2) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            return (int)a.size();
        }

        Interval* overlappingInterval(int pos, int pos2) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            Interval *res = new Interval;
            for (size_t i = 0; i < a.size(); ++i) {
                res->low = tree.start(a[i]);
                res->high = tree.end(a[i]);
                break;
            };
            return res;
        }

    private:
        IITree<int, int> tree;

};