#include <iostream>
#include <string>
#include <utility>
#include <deque>

#include "robin_map.h"
#include "robin_set.h"
#include "robin_hash.h"

#include "xxhash64.h"

#include "../htslib/sam.h"
#include "../htslib/hts.h"


int search_hts_file(char *infile , char *outfile, const int min_within_size, const int clip_length,
    const int threads) {
    // https://gist.github.com/PoisonAlien/350677acc03b2fbf98aa

    int result;

    samFile *fp_in = hts_open(infile, "r");
    result = hts_set_threads(fp_in, threads);
    if (result != 0) { exit(1); }

    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
    if (!bamHdr) {
        std::cerr << "Failed to read input header\nn";
        exit(1);
    }

    auto f_out = sam_open(outfile, "wbu");
    result = sam_hdr_write(f_out, bamHdr);
    if (result != 0) { exit(1); }

    const int max_scope = 100000;
    int total = 0;

    std::pair<uint64_t, bam1_t*> scope_item;
    std::deque<std::pair<uint64_t, bam1_t*>> scope;
    tsl::robin_set<uint64_t> read_names;

    while(sam_read1(fp_in,bamHdr,aln) > 0){

        if (scope.size() > max_scope) {
            scope_item = scope[0];

            if ( read_names.find(scope_item.first, scope_item.first) != read_names.end() ) {
                int result = sam_write1(f_out, bamHdr, scope_item.second);
                if (result == 0) { exit(1); }
                total += 1;
            }
            scope.pop_front();
        }

        uint16_t flag = aln->core.flag;

        if (flag & 1284 || aln->core.n_cigar == 0) {
            continue;
        }


//        std::string qname = std::string(bam_get_qname(aln), aln->core.l_qname);
//        const std::size_t precalculated_hash = std::hash<std::string>()(qname);

        const uint64_t precalculated_hash = XXHash64::hash(bam_get_qname(aln), aln->core.l_qname, 0);

        scope.push_back(std::make_pair(precalculated_hash, aln));

        if (read_names.find(precalculated_hash, precalculated_hash) == read_names.end()) {

            // Check for disfordant of supplementary
            if (~flag & 2 || flag & 2048) {
                read_names.insert(precalculated_hash);
                continue;
            }

            // Check for SA tag
            if (bam_aux_get(aln, "SA")) {
                read_names.insert(precalculated_hash);
                continue;
            }

            uint32_t *const cigar = bam_get_cigar(aln);

            // Check cigar
            for (uint32_t k = 0; k < aln->core.n_cigar; k++) {

                int op = bam_cigar_op(cigar[k]);
                int length = bam_cigar_oplen(cigar[k]);

                if ((op == BAM_CSOFT_CLIP ) && (length >= clip_length)) {  // || op == BAM_CHARD_CLIP
                    read_names.insert(precalculated_hash);
                    break;
                }
                if (( op == BAM_CINS || op == BAM_CDEL) && (length >= min_within_size)) {
                    read_names.insert(precalculated_hash);
                    break;
                }
            }
        }
	}

	while (scope.size() > 0) {
        scope_item = scope[0];
        if ( read_names.find(scope_item.first, scope_item.first) != read_names.end() ) {
            int result = sam_write1(f_out, bamHdr, scope_item.second);
            if (result == 0) { exit(1); }
            total += 1;
        }
        scope.pop_front();
    }

	bam_destroy1(aln);
	sam_close(fp_in);
	sam_close(f_out);

    return total;

};
