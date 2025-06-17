
#pragma once

#include <string>
#include <charconv>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <memory>
#include <utility>
#include <cstdint>
#include "superintervals.h"
#include "unordered_dense.h"


namespace Tr {

    template<typename Out>
    void split(const std::string &s, char delim, Out result) {
        size_t start = 0;
        size_t end = s.find(delim);
        while (end != std::string::npos) {
            if (start != end) {  // Skip empty strings inline
                *result++ = s.substr(start, end - start);
            }
            start = end + 1;
            end = s.find(delim, start);
        }
        if (start < s.length()) {  // Add the last segment
            *result++ = s.substr(start);
        }
    }

    std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        elems.reserve(std::count(s.begin(), s.end(), delim) + 1);  // Pre-allocate for expected size
        split(s, delim, std::back_inserter(elems));
        return elems;
    }

    bool endsWith(const std::string &mainStr, const std::string &toMatch) {
        if (mainStr.size() >= toMatch.size() &&
            mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
            return true;
        else
            return false;
    }

    enum FType {
         GFF3_NOI,
         GTF_NOI
    };

    class TrackBlock {
    public:
        std::string chrom, name, vartype, parent; // line
        int start{0}, end{0};
    };

    class GFFTrackBlock {
    public:
        std::string chrom, name, line, vartype;
        std::vector<std::string> parts;
        int start, end;
    };

    class Region {
    public:
        std::string chrom;
        int start, end;
    };

    /*
    * File reader. Non-indexed files are cached using TrackBlock items.
    */
    class TranscriptData {
    public:
        TranscriptData() = default;
        ~TranscriptData() = default;
        // The iterator state is cached here during iteration:
        std::string path;
        std::string chrom, chrom2, rid, vartype, parent;
        int start{-1}, stop{-1};
        int strand{0};
        int fileIndex;

        FType kind;

        std::shared_ptr<std::istream> fpu;
        std::string tp;

        int current_iter_index;
        int num_intervals;

        std::vector<std::string> parts;  // string split by delimiter
        std::vector<std::string> keyval;

        ankerl::unordered_dense::map< std::string, ankerl::unordered_dense::set<int64_t>> unique_gaps;
        ankerl::unordered_dense::map< std::string, SuperIntervals<int, std::pair<int, int>>> allBlocks;

        TrackBlock block;
        bool done{true};
        bool any_data{false};
		std::string variantString;
		std::vector<std::pair<int, int>> overlapping_tr_gaps;

        void open(const char* input_file) {
            std::string path = input_file;
            fileIndex = 0;
            done = false;

            if (endsWith(path, ".gff3")) {
                kind = GFF3_NOI;
            } else if (endsWith(path, ".gtf")) {
                kind = GTF_NOI;
            } else {
                std::cerr << "Error: only uncompressed GFF3 files are supported\n";
                std::cerr << "File name: " << path << std::endl;
                throw std::exception();
            }

            if (kind == GFF3_NOI || kind == GTF_NOI) {
                auto file_stream = std::make_shared<std::ifstream>(path);
                if (!file_stream->is_open()) {
                    std::cerr << "Error: opening track file " << path << std::endl;
                    throw std::runtime_error("Error opening file");
                }
                fpu = file_stream;
                std::string exon = "exon\t";

                ankerl::unordered_dense::map< std::string, std::vector<TrackBlock>> track_blocks;
                while (true) {
                    auto got_line = (bool)getline(*fpu, tp);
                    if (!got_line) {
                        done = true;
                        break;
                    }
                    fileIndex += 1;
                    if (tp[0] == '#') {
                        continue;
                    }
                    size_t col1_end = tp.find("\t");
                    size_t col2_end = tp.find("\t", col1_end + 1);
                    size_t col3_start = col2_end + 1;

                    // Keep only 'exon' lines
                    bool is_exon = tp.compare(col3_start, 5, exon) == 0;

                    if (!is_exon) {
                        continue;
                    }

                    TrackBlock b;
                    // Continue finding tab positions for other columns
                    size_t col3_end = tp.find("\t", col3_start);
                    size_t col4_start = col3_end + 1;
                    size_t col4_end = tp.find("\t", col4_start);
                    size_t col5_start = col4_end + 1;
                    size_t col5_end = tp.find("\t", col5_start);
                    size_t col6_start = col5_end + 1;
                    size_t col6_end = tp.find("\t", col6_start);
                    size_t col7_start = col6_end + 1;
                    size_t col7_end = tp.find("\t", col7_start);
                    size_t col8_start = col7_end + 1;
                    size_t col8_end = tp.find("\t", col8_start);
                    size_t col9_start = col8_end + 1;

                    // Extract needed fields directly
//                    b.line = tp;  // todo remove this
                    b.chrom = tp.substr(0, col1_end);

                    // Parse start and end using from_chars
                    int start_value = 0;
                    int end_value = 0;
                    std::from_chars(tp.data() + col4_start, tp.data() + col4_end, start_value);
                    std::from_chars(tp.data() + col5_start, tp.data() + col5_end, end_value);
                    b.start = start_value - 1;  // Use zero-based
                    b.end = end_value;

                    // Extract attributes (column 9)
                    std::string attributes = tp.substr(col9_start);
                    for (const auto &item : split(attributes, ';')) {
                        if (kind == GFF3_NOI) {
                            std::vector<std::string> keyval = split(item, '=');
                            if (keyval.size() >= 2) {  // Safety check
                                if (keyval[0] == "Parent") {
                                    b.parent = keyval[1];
                                    if (!b.name.empty()) {
                                        break;
                                    }
                                } else if (keyval[0] == "transcript_id") {
                                    b.name = keyval[1];
                                    if (!b.parent.empty()) {
                                        break;
                                    }
                                }
                            }
                        } else {  // GTF_NOI
                            // todo
                        }
                    }
                    if (b.name.empty()) {
                        continue;
                    }
                    track_blocks[b.name].push_back(std::move(b));
                }
                for (auto& kv : track_blocks) {
                    if (kv.second.size() <= 1) {
                        continue;
                    }
                    std::vector<Tr::TrackBlock>& blocks = kv.second;
                    std::sort(blocks.begin(), blocks.end(),
                        [](const TrackBlock& t1, const TrackBlock& t2) { return t1.start < t2.start; }
                    );
                    for (size_t i = 0; i < blocks.size() - 1; ++i) {
                        if (blocks[i].end < blocks[i+1].start) {
                            int64_t pos_key = (static_cast<int64_t>(blocks[i].end) << 32) | static_cast<int64_t>(blocks[i+1].start);
                            unique_gaps[blocks[i].chrom].insert(pos_key);
                        }
                    }
                }

            } else {
                std::cerr << "Error: file stype not supported for " << path << std::endl;
                throw std::exception();
            }

            for (const auto& item : unique_gaps) {
                for (auto& se : item.second ) {
                    int end = se & 0xFFFFFFFF;
                    int start = (se >> 32);
                    allBlocks[item.first].add(start, end, {start, end});
                }
                allBlocks[item.first].index();
                any_data = true;
            }
        }

        bool hasRefSkipGap(std::string& current_chrom, int start, int end, int tolerance=10) {
            overlapping_tr_gaps.clear();
            this->allBlocks[current_chrom].findOverlaps(start, end, this->overlapping_tr_gaps);
            for (const auto& ol : this->overlapping_tr_gaps) {
                if (std::abs(ol.first - start) < tolerance && std::abs(ol.second - end) < tolerance) {
                    return true;
                }
            }
            return false;
        }

        void writeUniqueGapsToBed(const char* out_file) {
            std::ofstream out(out_file);
            if (!out.is_open()) {
                throw std::runtime_error("Could not open file for writing");
            }
            for (const auto& item : unique_gaps) {
                for (auto& se : item.second ) {
                    int end = se & 0xFFFFFFFF;
                    int start = (se >> 32);
                    out << item.first << "\t" << start << "\t" << end << "\n";
                }
            }
            out.close();
        }
    };

} // namepsace Transc