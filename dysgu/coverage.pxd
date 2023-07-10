#cython: language_level=3
from dysgu.map_set_utils cimport CoverageTrack

cdef class GenomeScanner:
    """Takes care of scanning genome for reads and generating coverage track"""
    cdef public int mapq_threshold, max_cov, read_threads, buffer_size, clip_length, min_within_size, procs, buff_size, \
            current_cov, current_chrom, current_pos, reads_dropped, first, current_tid

    cdef public object input_bam, include_regions, regions_only, stdin, cov_track_path, overlap_regions, \
            staged_reads, current_bin, current_cov_array, depth_d, read_buffer, approx_read_length, last_tell, bam_iter
    cdef public bint paired_end, no_tell

    cdef CoverageTrack cpp_cov_track
