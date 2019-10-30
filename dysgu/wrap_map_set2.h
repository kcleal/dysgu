
#include <iostream>
#include <string>
//#include <tuple>
#include <utility>
//#include <functional>
//#include <unordered_map>
//#include <map>
#include "robin_map.h"
#include "robin_set.h"
#include "robin_hash.h"


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

            if(map.find(key) != map.end()) {
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


class Int2IntVecMap
{
    public:

        Int2IntVecMap() {}
        ~Int2IntVecMap() {}

        int size() { return map.size(); }

        void insert(int key, std::vector<int> value) { map[key] = value; }
        void erase(int key) { map.erase(key); }

        std::vector<int> get(int key) { return map.at(key); }

        int has_key(int key) {
            if (map.find(key) == map.end()) { return 0; }
            else { return 1; }
        }

    private:
        tsl::robin_map<int, std::vector<int>> map;

};

// https://stackoverflow.com/questions/10405030/c-unordered-map-fail-when-used-with-a-vector-as-key
// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
struct VectorHasher {
    std::size_t operator()(std::vector<int> const& vec) const {
      std::size_t seed = vec.size();
      for(auto& i : vec) {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
};


class IntVec2IntMap
{
    public:

        IntVec2IntMap() {}
        ~IntVec2IntMap() {}

        int size() { return map.size(); }

        void insert(std::vector<int> key, int value) { map[key] = value; }
        void erase(std::vector<int> key) { map.erase(key); }

        int get(std::vector<int> key) { return map.at(key); }

        int has_key(std::vector<int> key) {
            if (map.find(key) == map.end()) { return 0; }
            else { return 1; }
        }

    private:
        tsl::robin_map<std::vector<int>, int, VectorHasher> map;

};



