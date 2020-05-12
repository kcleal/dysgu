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

#include "./htslib/htslib/sam.h"
#include "./htslib/htslib/hfile.h"
#include "xxhash64.h"


int search_hts_alignments(char* infile, char* outfile, uint32_t min_within_size, uint32_t clip_length,
                      int threads) {

    int result;
    htsFile *fp_in = hts_open(infile, "r");
    if (threads > 1) {  // set additional threads
        result = hts_set_threads(fp_in, threads - 1);
        if (result != 0) { return -1; }
    }

    bam_hdr_t* samHdr = sam_hdr_read(fp_in);  // read header
    bam1_t* aln = bam_init1();  // initialize an alignment
    if (!samHdr) { return -1;}

    htsFile *f_out = hts_open(outfile, "wb0");
    result = sam_hdr_write(f_out, samHdr);
    if (result != 0) { return -1; }

    uint32_t max_scope = 100000;
    uint64_t total = 0;

    std::pair<uint64_t, bam1_t*> scope_item;
    std::deque<std::pair<uint64_t, bam1_t*>> scope;
    tsl::robin_set<uint64_t> read_names;

    uint16_t flag;
    uint32_t* cigar;
    bam1_t* bam_ptr;

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
        if (flag & 1284 || aln->core.n_cigar == 0 || aln->core.l_qname == 0) { continue; }

        const uint64_t precalculated_hash = XXHash64::hash(bam_get_qname(aln), aln->core.l_qname, 0);

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
            for (uint32_t k=0; k < aln->core.n_cigar; k++) {

                uint32_t op = bam_cigar_op(cigar[k]);
                uint32_t length = bam_cigar_oplen(cigar[k]);

                if ((op == BAM_CSOFT_CLIP ) && (length >= clip_length)) {  // || op == BAM_CHARD_CLIP
                    read_names.insert(precalculated_hash);
                    break;
                }

                if ((op == BAM_CINS || op == BAM_CDEL) && (length >= min_within_size)) {
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