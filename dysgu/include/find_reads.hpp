#pragma once

#include <cstdint>
#include <iostream>
#include <functional>
#include <string>
#include <utility>
#include <queue>
#include <map>
#include <fstream>
#include <sstream>

#include <htslib/sam.h>
#include <htslib/hfile.h>

#include "unordered_dense.h"
#include "xxhash64.h"
#include "transcripts.hpp"


void remove_extraneous_tags(bam1_t* src, int remove_extra_tags) {
    // needed NM, AS, SA, XS, HP, but remove some of the large offenders to save a bit of disk space
    int i=0;
    for (auto &t : {"ML", "MM", "MD", "RG", "f5", "st", "sn", "sv"}) {
        uint8_t* data = bam_aux_get(src, t);
        if (data != nullptr) {
            bam_aux_del(src, data);
        }
        i += 1;
        if (i >= remove_extra_tags) {
            break;
        }
    }
    uint8_t* qual = bam_get_qual(src);
    for (int i = 0; i < src->core.l_qseq; i++) {
        qual[i] = 60;
    }
}


class CoverageTrack
{
    public:

        CoverageTrack() {}
        ~CoverageTrack() {}

        std::vector<int32_t> cov_array;  // Assume coverage never overflows int32
        int max_coverage = 250;
        int index = 0;

        void set_max_cov(int m) {
            max_coverage = m;
        }

        void add(int index_start, int index_end) {
            cov_array[index_start / 10] += 1;
            cov_array[index_end / 10] -= 1;
        }

        int get_cov(int pos) {
            int idx = pos / 10;
            if (pos > 0 && idx <= cov_array.size()) {
                return cov_array[idx];
            }
            return -1;
        }

        bool cov_val_good(int current_tid, int aln_tid, int pos) {
            // Work out the coverage value up until alignment pos. Assumes that alignments dropping out of the queue
            // will not overlap the ones being read into the queue

            if (aln_tid != current_tid) { return true; }  // Assumed, but might not be the case

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
            for (int i = 0; i < cov_array.size(); i++) {

                if (i > index) {
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


int process_alignment(int& current_tid, std::string& current_chrom, std::deque<std::pair<uint64_t, bam1_t*>>& scope, std::vector<bam1_t*>& write_queue,
                       const int max_scope, const int max_write_queue, const int clip_length,
                       ankerl::unordered_dense::set<uint64_t>& read_names, CoverageTrack& cov_track, uint64_t& total,
                       const int n_chromosomes, const int mapq_thresh, const int check_clips, const int min_within_size,
                       sam_hdr_t **samHdr, htsFile **f_out,
                       std::string& temp_folder, const bool write_all, long *n_aligned_bases, int remove_extra_tags,
                       Tr::TranscriptData& transcript_gaps) {

    int result;
    bam1_t* aln;
    std::pair<uint64_t, bam1_t*> scope_item;

    if (scope.size() > max_scope) {
        scope_item = scope[0];
        aln = scope_item.second;
        // check if read is SV-read and in low coverage region, push to write queue
        if (read_names.find(scope_item.first) != read_names.end() && cov_track.cov_val_good(current_tid, aln->core.tid, aln->core.pos)) {
            write_queue.push_back(aln);
        } else {
            bam_destroy1(aln);
        }
        scope.pop_front();
    }

    // Check if write queue is full
    if (write_queue.size() > max_write_queue) {
        for (const auto& val: write_queue) {
            if (remove_extra_tags) {
                remove_extraneous_tags(val, remove_extra_tags);
            }
            result = sam_write1(*f_out, *samHdr, val);
            if (result < 0) { return -1; }
            total += 1;
            bam_destroy1(val);
        }
        write_queue.clear();
    }

    // Process current alignment as it arrives in scope. Add rname to hash if it is an SV read. Add coverage info
    aln = scope.back().second;
    const uint16_t flag = aln->core.flag;
    const uint32_t n_cigar = aln->core.n_cigar;

    // Skip uninteresting reads before putting on queue
    // unmapped, not primary, duplicate
    if (flag & 1284 || n_cigar == 0 || aln->core.l_qname == 0) {
        // Next item will overwrite this record
        return 0;
    }
    const uint16_t tid = aln->core.tid;
    if (tid != current_tid && tid <= n_chromosomes) {  // prepare coverage array
        const char* rname = sam_hdr_tid2name(*samHdr, current_tid);
        if (current_tid != -1 && rname != NULL) {  // Write last cov_track
            std::string out_name = temp_folder + "/" + std::string(rname) + ".dysgu_chrom.bin";
            cov_track.write_track(&out_name[0]);
        }
        rname = sam_hdr_tid2name(*samHdr, tid);
        current_chrom = rname;
        int chrom_length = sam_hdr_tid2len(*samHdr, tid);
        if (chrom_length == 0) {
            return 0;
        }
        cov_track.set_cov_array(chrom_length);
        current_tid = tid;
    }

    const uint32_t* cigar = bam_get_cigar(aln);

    if (aln->core.qual < mapq_thresh) {
        int index_start = aln->core.pos;
        for (uint32_t k=0; k < n_cigar; k++) {
            uint32_t op = bam_cigar_op(cigar[k]);
            uint32_t length = bam_cigar_oplen(cigar[k]);
            if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                index_start += length;
            } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                int index_end = index_start + length;
                cov_track.add(index_start, index_end);
                index_start = index_end;
                *n_aligned_bases += length;
            }
        }
        return 0;
    }

    size_t precalculated_hash = 0;
    precalculated_hash = XXHash64::hash(bam_get_qname(aln), aln->core.l_qname, precalculated_hash);
    if (!write_all) {
        precalculated_hash = XXHash64::hash(&aln->core.flag, sizeof(aln->core.flag), precalculated_hash);
    }

    scope.back().first = precalculated_hash;

    // Add a new item to the queue for next iteration
    scope.push_back(std::make_pair(0, bam_init1()));

    // add alignment to coverage track, check for indels and soft-clips
    bool sv_read = read_names.find(precalculated_hash) != read_names.end();
    if (sv_read) {
        return 0;
    }

    int index_start = aln->core.pos;
    for (uint32_t k=0; k < n_cigar; k++) {
        uint32_t op = bam_cigar_op(cigar[k]);
        uint32_t length = bam_cigar_oplen(cigar[k]);
        if (op == BAM_CDEL) {
            index_start += length;

        } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int index_end = index_start + length;
            cov_track.add(index_start, index_end);
            index_start = index_end;
            *n_aligned_bases += length;

        } else if (op == BAM_CREF_SKIP) {
            if (transcript_gaps.any_data) {
                int index_end = index_start + length;
                bool has_ref_skip_gap = transcript_gaps.hasRefSkipGap(current_chrom, index_start, index_end, 10);
                sv_read = (has_ref_skip_gap) ? false : true;
                if (has_ref_skip_gap) {
                    break;
                }
            }
            index_start += length;
        }
        if (!sv_read) {
            if ((check_clips) && (op == BAM_CSOFT_CLIP) && (length >= clip_length)) {
                sv_read = true;
            } else if ((op == BAM_CINS || op == BAM_CDEL) && (length >= min_within_size)) {
                sv_read = true;
            }
        }
    }

    if (!sv_read) { // not an sv read template yet
        // Check for discordant of supplementary
        if ((~flag & 2 && flag & 1) || flag & 2048) {
            sv_read = true;
        }
        // Check for SA tag
        else if (bam_aux_get(aln, "SA")) {
            sv_read = true;
        }
    }

    if (sv_read) {
        read_names.emplace(precalculated_hash);
    }
    return 0;
}


void collect_transcripts(const char* transcripts_file, const char* unique_gaps_file, Tr::TranscriptData& t_reader) {
    if (transcripts_file == nullptr || strlen(transcripts_file) == 0) {
        t_reader.done = true;
        return;
    }
    t_reader.open(transcripts_file);
    t_reader.writeUniqueGapsToBed(unique_gaps_file);
}

int search_hts_alignments(const char* infile, const char* outfile, uint32_t min_within_size, int clip_length, int mapq_thresh,
                          int threads, int paired_end, const char* temp_f, int max_coverage, const char* region,
                          const char* max_cov_ignore_regions, const char* fasta,
                          bool write_all, const char* write_mode,
                          const char* transcripts_file, const char* unique_gaps_file) {

    Tr::TranscriptData transcript_gaps = Tr::TranscriptData();
    collect_transcripts(transcripts_file, unique_gaps_file, transcript_gaps);

    const int check_clips = (clip_length > 0) ? 1 : 0;

    int result;

    samFile *fp_in = sam_open(infile, "r");
    result = hts_set_fai_filename(fp_in, fasta);
    if (result != 0) { return -1; }

    hts_idx_t *index;

    if (threads > 1) {  // set additional threads beyond main thread
        result = hts_set_threads(fp_in, threads - 1);
        if (result != 0) { return -1; }
    }

    sam_hdr_t* samHdr = sam_hdr_read(fp_in);  // read header

    int n_chromosomes = samHdr->n_targets;

    if (!samHdr) { return -1;}

    htsFile *f_out = hts_open(outfile, write_mode);  // "wb0"

    result = hts_set_threads(f_out, 1);
    if (result != 0) { return -1; }

    result = sam_hdr_write(f_out, samHdr);
    if (result != 0) { return -1; }

    // Reads need to be cached in order for coverage values to be accurately calculated on the fly
    int max_scope = 100000;
    int max_write_queue = 100000;
    if (paired_end == 0) {
        max_scope = 1000;
        max_write_queue = 1000;
    }
    uint64_t total = 0;

    std::pair<uint64_t, bam1_t*> scope_item;
    std::deque<std::pair<uint64_t, bam1_t*>> scope;
    std::vector<bam1_t*> write_queue;  // Write in blocks
    ankerl::unordered_dense::set<uint64_t> read_names;

    // Initialize first item in scope, set hash once read has been read
    scope.push_back(std::make_pair(0, bam_init1()));

    // Write coverage information using mosdepth style algorithm
    CoverageTrack cov_track;
    cov_track.max_coverage = max_coverage;

    std::string temp_folder = temp_f;

    int current_tid = -1;
    std::string current_chrom{};

    hts_itr_t *iter = NULL;

    std::string region_string = region;

    long n_aligned_bases = 0;
    int remove_extra_tags = 0;
    // iterate through regions
    if (region_string != ".,") {

        index = sam_index_load(fp_in,  infile);

        int string_pos;
        std::string delim = ",";
        std::string token;
        char *token_ptr;

        while ((string_pos = region_string.find(delim)) != std::string::npos) {

            token = region_string.substr(0, string_pos);
            token_ptr = &token[0];

            iter = sam_itr_querys(index, samHdr, token_ptr);
            if (iter == NULL) {
                std::cerr << "NULL iterator sam_itr_querys. Bad region format? " << region << std::endl;
                return -1;
            }
            region_string.erase(0, string_pos + delim.length());

            // Read alignment into the back of scope queue
            while ( sam_itr_next(fp_in, iter, scope.back().second) >= 0 ) {
                if (total == 0) {
                    for (auto &t : {"ML", "MM", "MD"}) {
                        uint8_t* data = bam_aux_get(scope.back().second, t);
                        if (data != nullptr) {
                            remove_extra_tags += 1;
                        } else {
                            break;
                        }
                    }
                }
                int success = process_alignment(current_tid, current_chrom, scope, write_queue,
                           max_scope, max_write_queue, clip_length,
                           read_names, cov_track, total,
                           n_chromosomes, mapq_thresh, check_clips, min_within_size,
                           &samHdr, &f_out, temp_folder, write_all, &n_aligned_bases, remove_extra_tags,
                           transcript_gaps);
                if (success < 0) {
                    std::cerr << "Failed to process input alignment. Stopping" << region << std::endl;
                    break;
                }
            }
        }
    } else {
        // iterate whole alignment file (no index needed)
        while (sam_read1(fp_in, samHdr, scope.back().second) >= 0) {
            if (total == 0) {
                for (auto &t : {"ML", "MM", "MD"}) {
                    uint8_t* data = bam_aux_get(scope.back().second, t);
                    if (data != nullptr) {
                        remove_extra_tags += 1;
                    } else {
                        break;
                    }
                }
            }
            int success = process_alignment(current_tid, current_chrom, scope, write_queue,
                       max_scope, max_write_queue, clip_length,
                       read_names, cov_track, total,
                       n_chromosomes, mapq_thresh, check_clips, min_within_size,
                       &samHdr, &f_out, temp_folder, write_all, &n_aligned_bases, remove_extra_tags,
                       transcript_gaps);
            if (success < 0) {
                std::cerr << "Failed to process input alignment. Stopping" << region << std::endl;
                break;
            }
        }
    }

    // Deal with reads still in the queues
    bam1_t* aln;
    while (scope.size() > 0) {
        scope_item = scope[0];
        aln = scope_item.second;
        if (read_names.find(scope_item.first) != read_names.end() && cov_track.cov_val_good(current_tid, aln->core.tid, aln->core.pos)) {
            write_queue.push_back(aln);
        } else {
            bam_destroy1(aln);
        }
        scope.pop_front();
    }

    for (const auto& val: write_queue) {
        if (remove_extra_tags) {
            remove_extraneous_tags(val, remove_extra_tags);
        }
        result = sam_write1(f_out, samHdr, val);
        if (result < 0) { return -1; }
        total += 1;
        bam_destroy1(val);
    }

    // write last chrom to coverage track
    if (current_tid >= 0 && current_tid <= n_chromosomes) { // tid can sometimes be 65535
        const char* rname = sam_hdr_tid2name(samHdr, current_tid);
        if (rname != NULL) {
            std::string out_name = temp_folder + "/" + std::string(rname) + ".dysgu_chrom.bin";
            cov_track.write_track(&out_name[0]);
        }
    }

    std::ofstream bases_file;
    bases_file.open(temp_folder + "/n_aligned_bases.txt");
    bases_file << n_aligned_bases << std::endl;
    bases_file.close();

    result = hts_close(fp_in);
    if (result != 0) { return -1; };

    result = hts_close(f_out);
    if (result < 0) { return -1; };

    f_out = NULL;

    return total;

}
