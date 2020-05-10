#include <cstdint>
#include <iostream>
#include <functional>
#include <string>
#include <utility>
#include <queue>
#include <map>
#include "robin_map.h"
#include "robin_set.h"
#include "robin_hash.h"
#include "../htslib/htslib/sam.h"
#include "../htslib/htslib/hfile.h"
#include "xxhash64.h"

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


int search_hts_alignments(char* infile, char* outfile, uint32_t min_within_size, uint32_t clip_length,
                      int threads) {

    int result;
    htsFile *fp_in = hts_open(infile, "r");
    if (threads > 0) {
        result = hts_set_threads(fp_in, threads);
        if (result != 0) { return -1; }
    }

    bam_hdr_t* samHdr = sam_hdr_read(fp_in);  // read header
    bam1_t* aln = bam_init1();  // initialize an alignment
    if (!samHdr) { return -1;}

    htsFile *f_out = hts_open(outfile, "wb0");
    result = sam_hdr_write(f_out, samHdr);
    if (result != 0) { return -1; }

    XXHash64::XXHash64 hasher(0);

    uint32_t max_scope = 100000;
    uint64_t total = 0;

    std::pair<uint64_t, bam1_t*> scope_item;
    std::deque<std::pair<uint64_t, bam1_t*>> scope;
    tsl::robin_set<uint64_t> read_names;

    uint16_t flag;
    uint32_t* cigar;
//    uint32_t k, op, length;
//    uint64_t precalculated_hash;
    bam1_t* bam_ptr;
//    bam1_t bamt;

    while (sam_read1(fp_in, samHdr, aln) >= 0) {

        if (scope.size() > max_scope) {
            scope_item = scope[0];

            if (read_names.find(scope_item.first, scope_item.first) != read_names.end()) {
                result = sam_write1(f_out, samHdr, scope_item.second);
                if (result < 0) { return -1; }
                total += 1;
            }
            // free(dereference(scope_item.second).data)
            // free(scope_item.second)
            bam_destroy1(scope_item.second);
            scope.pop_front();
        }

        // Process current alignment
        flag = aln->core.flag;
        if (flag & 1284 or aln->core.n_cigar == 0 || aln->core.l_qname == 0) { continue; }

//        precalculated_hash = XXHash64::hash(bam_get_qname(aln), &aln->core.l_qname, 0);
        hasher.add(bam_get_qname(aln), aln->core.l_qname); // call add() as often as you like to ...
        // and compute hash:
        uint64_t precalculated_hash = hasher.hash();

        bam_ptr = bam_dup1(aln);
        scope.push_back(std::make_pair(precalculated_hash, bam_ptr));  // bam_dup1(aln)

        if (read_names.find(precalculated_hash, precalculated_hash) == read_names.end()) {

            // Check for discordant of supplementary
            if (~flag & 2 || flag & 2048) {
                read_names.insert(precalculated_hash);
                continue;
            }
            // Check for SA tag
            if (bam_aux_get(aln, "SA")) {
                read_names.insert(precalculated_hash);
                continue;
            }

            cigar = bam_get_cigar(aln);

            // Check cigar
            for (uint32_t k = 0; k < aln->core.n_cigar; k++) {
//            for k in range(0, dereference(aln).core.n_cigar):

                uint32_t op = bam_cigar_op(cigar[k]);
                //if (!op) { break; }

                uint32_t length = bam_cigar_oplen(cigar[k]);
                //if (!length) { break; }

                if ((op == BAM_CSOFT_CLIP ) && (length >= clip_length)) {  // || op == BAM_CHARD_CLIP
                    read_names.insert(precalculated_hash);
                    break;
                }

                if ((op == BAM_CINS or op == BAM_CDEL) && (length >= min_within_size)) {
                    read_names.insert(precalculated_hash);
                    break;
                }
            }
        }
    }

    while (scope.size() > 0) {
        scope_item = scope[0];
        if (read_names.find(scope_item.first, scope_item.first) != read_names.end()) {
            result = sam_write1(f_out, samHdr, scope_item.second);
            if (result < 0) { return -1; };
            total += 1;
        }

        scope.pop_front();
    }

    result = hts_close(fp_in);
    if (result != 0) { return -1; };

    result = hts_close(f_out);
    if (result < 0) { return -1; };

    f_out = NULL;

    bam_destroy1(aln);
    return total;

}

