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
//#include "htslib/sam.h"
//#include "htslib/hfile.h"
#include "xxhash64.h"


class CoverageTrack
{
    public:

        CoverageTrack() {}
        ~CoverageTrack() {}

        std::vector<int32_t> cov_array;  // Assume coverage never overflows int32
        int max_coverage;
        int index = 0;

        void add(int index_start, int index_end) {
            cov_array[index_start / 10] += 1;
            cov_array[index_end / 10] -= 1;
        }

        bool cov_val_good(bam1_t* aln, int current_tid) {
            // Work out the coverage value up until alignment pos. Assumes that alignments dropping out of the queue
            // will not overlap the ones being read into the queue

            if (aln->core.tid != current_tid) {
                return true;  // Assumed, but might not be the case
            }

            int pos = aln->core.pos;
            int target_index = (pos / 10) + 1;  // Make sure the pos bin is tested
            int current_cov = cov_array[index];

            if (target_index > index + 1) {
                // Cumulative sum of the cov array up until start of the alignment
                for (int i=index + 1; i < target_index; i++) {
                    current_cov += cov_array[i];
                    cov_array[i] = current_cov;
                    index += 1;
                }
            }

            assert(current_cov >= 0);

            if (current_cov > max_coverage) {
                return false;
            }
//            std::cerr << current_cov << std::endl;
//            exit(0);
            return true;
        }

        void set_cov_array(int chrom_length) {
            cov_array.clear();
            cov_array.resize((chrom_length / 10) + 1, 0);
            index = 0;
        }

        void write_track(char* out_name) {

            std::ofstream file_out (out_name, std::ios::binary);

            int current_cov = 0;
            int16_t high = 32000;

            //for (const auto& v: cov_array) {
            for (int i = 0; i < cov_array.size(); i++) {

                if (i >= index) {
                    int v = cov_array[i];
                    current_cov += v;
                } else {
                    current_cov = cov_array[i];
                }

                if (current_cov > high) {
                    file_out.write((char *)&high, sizeof(int16_t));
                } else {
                    file_out.write((char *)&current_cov, sizeof(int16_t));
                }
            }

            file_out.close();
        }
};


int search_hts_alignments(char* infile, char* outfile, uint32_t min_within_size, int clip_length, int mapq_thresh,
                          int threads, int paired_end, char* temp_f, int max_coverage) {

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
     cov_track.max_coverage = max_coverage;

     std::string temp_folder = temp_f;

    int current_tid = -1;
    bam1_t* aln;

    // Read alignment into the back of scope queue
    while (sam_read1(fp_in, samHdr, scope.back().second) >= 0) {

        if (scope.size() > max_scope) {
            scope_item = scope[0];

            // check if read is SV-read and in low coverage region, push to write queue
            if (read_names.find(scope_item.first) != read_names.end() && cov_track.cov_val_good(scope_item.second, current_tid)) {
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
        aln = scope.back().second;
        const uint16_t flag = aln->core.flag;

        // Skip uninteresting reads before putting on queue
        // unmapped, not primary, duplicate
        if (flag & 1284 || aln->core.n_cigar == 0 || aln->core.l_qname == 0 || aln -> core.qual < mapq_thresh ) {
            // Next item will overwrite this record
            continue;
        }

        const uint16_t tid = aln->core.tid;

        if (tid != current_tid && tid <= n_chromosomes) {  // prepare coverage array
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
//        bool to_check = true;
//        if (aln -> core.qual < mapq_thresh) {
//            to_check = false;
//        }

        if (read_names.find(precalculated_hash) != read_names.end()) { sv_read = true; }

        int index_start = aln->core.pos;

        for (uint32_t k=0; k < aln->core.n_cigar; k++) {

            uint32_t op = bam_cigar_op(cigar[k]);
            uint32_t length = bam_cigar_oplen(cigar[k]);

            if (!sv_read) {
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

        if (!sv_read) { //&& (read_names.find(precalculated_hash) == read_names.end())) { // not an sv read template yet

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
        if (read_names.find(scope_item.first) != read_names.end() && cov_track.cov_val_good(scope_item.second, current_tid)) {
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