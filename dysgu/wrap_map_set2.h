#include <cstdint>
#include <iostream>
#include <functional>
#include <string>
#include <utility>
#include <queue>
#include <map>
#include <cassert>
#include "robin_map.h"
#include "robin_set.h"
#include "robin_hash.h"

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

            // Look for v
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

            for (auto& val: outList[u]) {
                if (val.first == v) {
                    val.second += w;
                    return;
                }
            }

            outList[u].push_back(std::make_pair(v, w));
            inList[v].push_back(std::make_pair(u, w));
            n_edges += 1;
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


        std::vector<int> neighbors(int u) {  // Get all out edges

            std::vector<int> neigh;

            for (const auto& val: outList[u]) {
                neigh.push_back(val.first);
            }

            return neigh;

        }


        std::vector<PairW2> forInEdgesOf(int u) {

            std::vector<PairW2> inEdges;

            for (const auto& val: inList[u]) {
                inEdges.push_back(val);
            }

            return inEdges;

        }


};


class SimpleGraph {
    // Undirected, weighted
    public:

        SimpleGraph() {}  // Construct
        ~SimpleGraph() {}

        // Used for undirected graph
        std::vector<std::vector<PairW>> adjList;

        // Used for directed graph
//        std::vector<std::vector<PairW>> inList;
//        std::vector<std::vector<PairW>> outList;

        int N = 0;
        int n_edges = 0;


        int addNode() {

            int n = N;

            std::vector<PairW> node;
//            node.push_back(std::make_pair(-1, -1));
            adjList.push_back(node);

            N++;
            return n;
        }

        int hasEdge(int u, int v) {
            if ((u > N ) || (v > N)) { return 0; }

            // Look for v
            for (const auto& val: adjList[u]) {
                if (val.first == v) {
                    return 1;
                }
            }
            return 0;

        }

        void addEdge(int u, int v, uint8_t w) {

            // Assume hasEdge has been call
            // If node has no edge add first in, otherwise push a new node
//            if ((adjList[u].size() == 1) && (adjList[u].front().first == -1)) {
//                adjList[u].front().first = v;
//                adjList[u].front().second = w;
//            } else {
//                adjList[u].push_back(std::make_pair(v, w));
//            };
//
//            if ((adjList[v].size() == 1) && (adjList[v].front().first == -1)) {
//                adjList[v].front().first = u;
//                adjList[v].front().second = w;
//            } else {
//                adjList[v].push_back(std::make_pair(u, w));
//            };

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

        std::vector<int> neighbors(int u) {

            std::vector<int> neigh;

            for (const auto& val: adjList[u]) {
                neigh.push_back(val.first);
            }

            return neigh;

        }

        void removeNode(int u) {
            // Set node edges to empty
            std::vector<PairW> node;
            node.push_back(std::make_pair(-1, -1));
            adjList[u] = node;
        }

        std::vector<int> connectedComponents() {

            std::vector<bool> visited(adjList.size(), false);

            std::vector<int> components;

            for (int u=0; u<N; u++) {

                if (visited[u] == false) {


                    components.push_back(u);
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
                            components.push_back(v);
                            visited[v] = true;

                            for (const auto& val2: adjList[v]) {

                                if ((val2.first != -1) && (visited[val2.first] == false)) {
                                    queue.push_back(val2.first);
                                }
                            }
                        }
                    }

                    components.push_back(-1);  // -1 is end of component
                }
            }

            return components;


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
        tsl::robin_map<int, int> map;

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
        tsl::robin_set<int> set;

};


class StrSet
{
    public:

        StrSet() {}
        ~StrSet() {}

        int size() { return set.size(); }

        void insert(std::string key) { set.insert(key); }
        void erase(std::string key) { set.erase(key); }

        int has_key(std::string key) {
            if (set.find(key) == set.end()) { return 0; }
            else { return 1; }
        }

    private:
        tsl::robin_set<std::string> set;

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

        void add_tuple_key(uint64_t packed_data, int index) {
            assert (index == index_key.size());
            string_key[packed_data] = index;
            index_key_map.push_back(packed_data);
        }
           int has_tuple_key(uint64_t packed_data) {

            tsl::robin_map<uint64_t, int>::const_iterator got = string_key.find(packed_data, packed_data);
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

        std::vector<int> idx_2_vec(int index) {
            std::vector<int> v = {0, 0, 0, 0};
            uint64_t packed_data = index_key_map[index];
            v[0] = packed_data & 15;  // 4 bits set to 1
            v[1] = (packed_data >> 4) & 4294967295;  // 1 x 32 bits
            v[2] = (packed_data >> 36) & 8388607;  // 1 x 23 bits
            v[3] = (packed_data >> 59) & 15;

            return v;
        }

    private:

        uint64_t last_key = 0;  // Set this on call to has_tuple_key to prevent calculating twice
        int last_index = 0;

        tsl::robin_map<uint64_t, int> string_key;
        std::vector<uint64_t> index_key_map;

};