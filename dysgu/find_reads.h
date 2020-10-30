#include <cstdint>
#include <iostream>
#include <functional>
#include <string>
#include <utility>
#include <queue>
#include <map>
#include <fstream>

#include "robin_hood.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/hfile.h"
#include "xxhash64.h"


class CoverageTrack
{
    public:

        CoverageTrack() {}
        ~CoverageTrack() {}

        std::vector<int32_t> cov_array;

        void add(int index_start, int index_end) {
            int i = index_start / 10;
//            if (index_start >= 0 && index_end < cov_array.size())
            //if (cov_array[i] < 32000 && cov_array[i] > -32000) {
            cov_array[i] += 1;
            cov_array[index_end / 10] -= 1;
            //}
        }

        void set_cov_array(int chrom_length) {
            cov_array.clear();
            cov_array.resize((chrom_length / 10) + 1, 0);
        }

        void write_track(char* out_name) {

            std::ofstream file_out (out_name, std::ios::binary);

            int current_cov = 0;
            int16_t high = 32000;
//            int n_over_high = 0;
            for (const auto& v: cov_array) {
//                if (v != 0) {
                current_cov += v;
//                }
                if (current_cov > high) {
                    file_out.write((char *)&high, sizeof(int16_t));
//                    n_over_high += 1;
                } else {
                    file_out.write((char *)&current_cov, sizeof(int16_t));
                }
            }
//            if (n_over_high > 0) {
//                std::cerr << "N coverage bins out of range: " << n_over_high << std::endl;
//            }
            file_out.close();
        }
};


int search_hts_alignments(char* infile, char* outfile, uint32_t min_within_size, int clip_length, int mapq_thresh,
                          int threads, int paired_end, char* temp_f) {

    const int check_clips = (clip_length > 0) ? 1 : 0;

    int result;
    htsFile *fp_in = hts_open(infile, "r");

    if (threads > 1) {  // set additional threads beyond main thread
        result = hts_set_threads(fp_in, threads - 1);
        if (result != 0) { return -1; }
    }

    bam_hdr_t* samHdr = sam_hdr_read(fp_in);  // read header

    int n_chromosomes = samHdr->n_targets;

    if (!samHdr) { return -1;}

    htsFile *f_out = hts_open(outfile, "wb0");
    result = hts_set_threads(f_out, 1);
        if (result != 0) { return -1; }

    result = sam_hdr_write(f_out, samHdr);
    if (result != 0) { return -1; }

    int max_scope = 100000;
    int max_write_queue = 100000;

    if (paired_end == 0) {
        max_scope = 100;
        max_write_queue = 100;
    }
    uint64_t total = 0;

    std::pair<uint64_t, bam1_t*> scope_item;
    std::deque<std::pair<uint64_t, bam1_t*>> scope;
    std::vector<bam1_t*> write_queue;  // Write in blocks
    robin_hood::unordered_set<uint64_t> read_names;

    // Initialize first item in scope, set hash once read has been read
    scope.push_back(std::make_pair(0, bam_init1()));

    // Write coverage information using mosdepth style algorithm
     CoverageTrack cov_track;
     std::string temp_folder = temp_f;

    int current_tid = -1;

    // Read alignment into the back of scope queue
    while (sam_read1(fp_in, samHdr, scope.back().second) >= 0) {

        if (scope.size() > max_scope) {
            scope_item = scope[0];

            if (read_names.find(scope_item.first) != read_names.end()) {
                write_queue.push_back(scope_item.second);
            } else {
                bam_destroy1(scope_item.second);
            }
            scope.pop_front();
        }

        // Check if write queue is full
        if (write_queue.size() > max_write_queue) {
            for (const auto& val: write_queue) {
                result = sam_write1(f_out, samHdr, val);
                if (result < 0) { return -1; }
                total += 1;
                bam_destroy1(val);
            }
            write_queue.clear();
        }

        // Process current alignment as it arrives in scope. Add rname to hash if it is an SV read. Add coverage info
        bam1_t* aln = scope.back().second;
        const uint16_t flag = aln->core.flag;

        // Skip uninteresting reads before putting on queue
        // unmapped, not primary, duplicate
        if (flag & 1284 || aln->core.n_cigar == 0 || aln->core.l_qname == 0 || aln -> core.qual < mapq_thresh) {
            // Next item will overwrite this record
            continue;
        }

        const uint16_t tid = aln->core.tid;

        if (tid != current_tid & tid <= n_chromosomes) {  // prepare coverage array
            if (current_tid != -1 ) {
                const char* rname = sam_hdr_tid2name(samHdr, current_tid);
                if (rname != NULL) {
//                    std::cerr << rname << std::endl;
                    std::string out_name = temp_folder + "/" + std::string(rname) + ".dysgu_chrom.bin";
                    cov_track.write_track(&out_name[0]);
                }
            }

            int chrom_length = sam_hdr_tid2len(samHdr, tid);
            if (chrom_length == 0) {
                continue;
            }
            cov_track.set_cov_array(chrom_length);
            current_tid = tid;
        }

        const uint64_t precalculated_hash = XXHash64::hash(bam_get_qname(aln), aln->core.l_qname, 0);
        scope.back().first = precalculated_hash;

        // Add a new item to the queue for next iteration
        scope.push_back(std::make_pair(0, bam_init1()));

        // add alignment to coverage track, check for indels and soft-clips
        const uint32_t* cigar = bam_get_cigar(aln);
        bool sv_read = false;

        if (read_names.find(precalculated_hash) != read_names.end()) { sv_read = true; }

        int index_start = aln->core.pos;

        for (uint32_t k=0; k < aln->core.n_cigar; k++) {

            uint32_t op = bam_cigar_op(cigar[k]);
            uint32_t length = bam_cigar_oplen(cigar[k]);

            if (~sv_read) {
                if ((check_clips) && (op == BAM_CSOFT_CLIP) && (length >= clip_length)) {  // || op == BAM_CHARD_CLIP
                    read_names.insert(precalculated_hash);
                    sv_read = true;

                } else if ((op == BAM_CINS || op == BAM_CDEL) && (length >= min_within_size)) {
                    read_names.insert(precalculated_hash);
                    sv_read = true;
                }
            }

            if (op == BAM_CDEL) {
                index_start += length;

            } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                int index_end = index_start + length;
                cov_track.add(index_start, index_end);
                index_start = index_end;
            }
        }

        if (~sv_read) { //&& (read_names.find(precalculated_hash) == read_names.end())) { // not an sv read template yet

            // Check for discordant of supplementary
            if ((~flag & 2 && flag & 1) || flag & 2048) {
                read_names.insert(precalculated_hash);
                continue;
            }
            // Check for SA tag
            if (bam_aux_get(aln, "SA")) {
                read_names.insert(precalculated_hash);
                continue;
            }
        }
    }

    while (scope.size() > 0) {
        scope_item = scope[0];
        if (read_names.find(scope_item.first) != read_names.end()) {
            write_queue.push_back(scope_item.second);
        } else {
            bam_destroy1(scope_item.second);
        }
        scope.pop_front();
    }

    for (const auto& val: write_queue) {
        result = sam_write1(f_out, samHdr, val);
        if (result < 0) { return -1; }
        total += 1;
        bam_destroy1(val);
    }

    // write last chrom to coverage track

    if (current_tid >= 0 && current_tid <= n_chromosomes) { // tid can sometimes be 65535
        const char* rname = sam_hdr_tid2name(samHdr, current_tid);
        if (rname != NULL) {
//            std::cerr << rname << std::endl;
            std::string out_name = temp_folder + "/" + std::string(rname) + ".dysgu_chrom.bin";
            cov_track.write_track(&out_name[0]);
        }
    }

    result = hts_close(fp_in);
    if (result != 0) { return -1; };

    result = hts_close(f_out);
    if (result < 0) { return -1; };

    f_out = NULL;

    return total;

}