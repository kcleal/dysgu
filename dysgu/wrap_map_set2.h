
#include <iostream>
#include <functional>
#include <string>
#include <utility>
#include <queue>
#include <map>
#include "robin_map.h"
#include "robin_set.h"
#include "robin_hash.h"





class DiGraph {
    // Directed non-weighted
    public:

        DiGraph() {}  // Construct
        ~DiGraph() {}

        std::vector<std::vector<int>> inList;
        std::vector<std::vector<int>> outList;

        int N = 0;
        int n_edges = 0;

        int addNode() {

            int n = N;

            std::vector<int> node;
            std::vector<int> in_node;

            outList.push_back(node);  // empty, no edges
            inList.push_back(in_node);
            N++;
            return n;
        }

        int hasEdge(int u, int v) {
            if ((u > N ) || (v > N)) { return 0; }

            // Look for v
            for (const auto& val: outList[u]) {
                if (val == v) {
                    return 1;
                }
            }
            return 0;

        }

        void addEdge(int u, int v) {
            // Assume hasEdge has been call
            outList[u].push_back(v);
            inList[v].push_back(u);
            n_edges += 1;
        }


        int numberOfNodes() { return N; }


        std::vector<int> neighbors(int u) {  // Get all out edges

            std::vector<int> neigh;

            for (const auto& val: outList[u]) {
                neigh.push_back(val);
            }

            return neigh;

        }


        std::vector<int> forInEdgesOf(int u) {

            std::vector<int> inEdges;

            for (const auto& val: inList[u]) {
                inEdges.push_back(val);
            }

            return inEdges;

        }


};


typedef std::pair<int, int> PairW;  // v, weight


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

        void addEdge(int u, int v, int w) {

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




//class PairScope {
//
//    public:
//
//        int max_dist;
//        int n_references;
//        int local_chrom = -1;
//
//        std::vector<std::map<int, std::pair<int, int>>> forward_scope;
//        std::vector<int> in_scope;
//        std::deque<std::pair<int, int>> local_scope;
//
//        PairScope() {}
//        ~PairScope() {}
//
//        void add_params(int m, int n) {
//            max_dist = m;
//            n_references = n;
//            for (int i = 0; i < n_references + 1; i++) {
//                std::map<int, std::pair<int, int>> scope_item;
//                forward_scope.push_back(scope_item);
//                in_scope.push_back(0);
//            }
//        }
//
//        void add_item(int c, int p, int node_name, int c_pos) {
//
//            forward_scope.at(c)[p] = std::make_pair(node_name, c_pos);
//        }
//
//        void drop_one(int c) {
//            in_scope[c] -= 1;
//            if (in_scope[c] == 0) {
//                forward_scope[c].clear();
//            }
//        }
//
//        std::vector<int> update(int node_name, int c_chrom, int c_pos, int chrom2, int pos2) {
//
//            std::vector<int> out_v;
//            std::pair<int, int> int_pair;
//
//            if (c_chrom != local_chrom) {
//                local_chrom = c_chrom;
//                local_scope.clear();
//
//                local_scope.push_back(std::make_pair(c_pos, chrom2));
//                add_item(chrom2, pos2, node_name, c_pos);
//
//                return out_v;
//            }
//
//            // Clear out of scope items
//            while (local_scope.size() > 0) {
//                int_pair = local_scope[0];
//
//                if (std::abs(int_pair.first - c_pos) > max_dist) {
//                    local_scope.pop_front();
//                    drop_one(int_pair.second);
//                } else {
//                    break;
//                }
//
//            }
//
//            if (forward_scope.at(chrom2).size() == 0) {
//                forward_scope.at(chrom2)[pos2] = std::make_pair(node_name, c_pos);
//                return out_v;
//            }
//
//            // Find closest one or two keys in map
//            // https://stackoverflow.com/questions/28404971/finding-the-closest-or-exact-key-in-a-stdmap
//
//            std::map<int, std::pair<int, int>> map;
//            map = forward_scope.at(chrom2);
//            std::map<int, std::pair<int, int>>::iterator low, prev;
//
//            low = map.lower_bound(pos2);
//
//            int item_forward_pos;
//            int item_node;
//            int item_local_pos;
//
//
//            if (low == map.end()) {
//
//                prev = std::prev(low);
//                item_forward_pos = prev->first;
//                item_node = prev->second.first;
//                item_local_pos = prev->second.second;
////                std::cerr << "here" << item_node << std::endl;
////                std::cerr <<  pos2 << " " << node_name << std::endl;
////                std::cerr <<  item_forward_pos << " " << item_node << std::endl;
//
////                for(std::map<int, std::pair<int, int>>::iterator it = map.begin(); it != map.end(); ++it) {
////
////                    std::cerr << it->first << "\n";
////                }
//
//                if ((std::abs(item_local_pos - c_pos) < max_dist) and (std::abs(item_forward_pos - pos2) < max_dist)) {
//                    out_v.push_back(item_node);
//
//                }
//
//            } else if (low == map.begin()) {
//                // First element
//                item_forward_pos = low->first;
//                item_node = low->second.first;
//                item_local_pos = low->second.second;
//
//                if (std::abs(item_local_pos - c_pos) < max_dist && std::abs(item_forward_pos - pos2) < max_dist) {
//                    out_v.push_back(item_node);
//                }
//
//            } else {
//                // In middle of map, add lower and higher elements
//
//                prev = std::prev(low);
//                item_forward_pos = prev->first;
//                item_node = prev->second.first;
//                item_local_pos = prev->second.second;
//
//                if (std::abs(item_local_pos - c_pos) < max_dist && std::abs(item_forward_pos - pos2) < max_dist) {
//                    out_v.push_back(item_node);
//                }
//
//                item_forward_pos = low->first;
//                item_node = low->second.first;
//                item_local_pos = low->second.second;
//
//                if (std::abs(item_local_pos - c_pos) < max_dist && std::abs(item_forward_pos - pos2) < max_dist) {
//                    out_v.push_back(item_node);
//                }
//
//            }
//
//            local_scope.push_back(std::make_pair(c_pos, chrom2));
//            forward_scope.at(chrom2)[pos2] = std::make_pair(node_name, c_pos);
//
//            return out_v;
//         }
//
//
//};











