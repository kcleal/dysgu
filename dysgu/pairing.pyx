#cython: language_level=3, boundscheck=False, wraparound=False
#distutils: language=c++

"""
Find the best path through a collection of alignments
Works with paried-end reads or single contigs. Borrows the pairing heuristic from bwa.
"""
from __future__ import absolute_import

import click
import array
from cpython cimport array
import numpy as np
cimport numpy as np
from libc.math cimport exp, log, sqrt, fabs
from libcpp.vector cimport vector


#DTYPE = np.int
#ctypedef np.int_t DTYPE_t
DTYPE = np.float
ctypedef np.float_t DTYPE_t


cdef float erfcc(float x):
    """Complementary error function."""
    cdef float z, t, r, p1, p2, p3, p4, p5, p6, p7, p8, p9
    z = fabs(x)
    t = (1. / (1. + 0.5*z))
    p1 = -.82215223+t*.17087277
    p2 = 1.48851587+t*p1
    p3 = -1.13520398+t*p2
    p4 = .27886807+t*p3
    p5 = -.18628806+t*p4
    p6 = .09678418+t*p5
    p7 = .37409196+t*p6
    p8 = 1.00002368+t*p7
    p9 = 1.26551223+t*p8
    r = t * exp(-z*z-p9)
    if x >= 0.:
        return r
    else:
        return 2. - r


cdef float normcdf(float x, float mu, float sigma):
    cdef float t, y
    t = x-mu
    y = 0.5*erfcc(-t/(sigma*sqrt(2.0)))
    if y > 1.0:
        y = 1.0
    return y


cdef int is_proper_pair(float d, float pos1, float pos2, float strand1, float strand2, float mu, float sigma):

    if d <= (mu + 4*sigma):
        if pos1 < pos2 and strand1 == 1. and strand2 == -1.:
            return 1
        if pos2 < pos1 and strand2 == 1. and strand1 == -1.:
            return 1
    return 0


cdef float bwa_pair_score(float d, float mu, float sigma, float match_score, float u):
    """
    Calculate the bwa pairing cost
    :param d: seperation distance
    :param mu: insert size mean
    :param sigma: insert size std
    :param match_score: score gained from a matched base
    :param u: constant parameter, worst cost possible
    :return: The bwa pairing cost (float), whether the read is FR (int)
    """
    cdef float prob, c
    prob = (1 - normcdf(d, mu, sigma))
    c = -match_score * (log(prob)/log(4))
    if c < u:
        return c
    else:
        return u

cdef tuple optimal_path(np.ndarray[DTYPE_t, ndim=2] segments,
                        float contig_length,
                        float mu,
                        float sigma,
                        float max_insertion,
                        float min_aln,
                        float max_homology,
                        float ins_cost,
                        float hom_cost,
                        float inter_cost,
                        float U,
                        float match_score,
                        ):
    """
    The scoring has become quite complicated.
     = Last alignment score - jump cost

    where each jump has a cost:
    jump cost = S + microhomology_length + insertion_length        # S = 10 for intra, 20 for inter

    The last score takes into account mismatched and gaps

    :param mapped:  The input dataframe with candidate alignments
    :param contig_length:   Length of the contig
    :param max_insertion: Arbitrary high number
    :param min_aln: The minimum sequence which is not part of a microhomology overlap
    :param max_homology: Arbitrary high number
    :return: tuple
    """

    # Start at first node then start another loop running backwards through the preceeding nodes.
    # Choose the best score out of the in edges.
    # Use a special start and end node to score the first and last alignments

    cdef np.ndarray[np.int_t, ndim=1] pred = np.zeros(segments.shape[0], dtype=np.int)
    cdef np.ndarray[np.float_t, ndim=1] node_scores = np.zeros(segments.shape[0], dtype=np.float)
    # Next best node score, for finding the secondary path
    cdef np.ndarray[np.float_t, ndim=1] nb_node_scores = np.zeros(segments.shape[0], dtype=np.float)

    normal_jumps = set([])  # Keep track of which alignments form 'normal' pairs between read-pairs. Alas native python

    cdef int i, j, p, normal_end_index, proper_pair # FR,
    cdef float chr1, pos1, start1, end1, score1, row_index1, strand1, r1,\
               chr2, pos2, start2, end2, score2, row_index2, strand2, r2, \
               micro_h, ins, best_score, next_best_score, best_normal_orientation, current_score, total_cost,\
               S, sc, max_s, path_score, jump_cost, distance

    # Deal with first score
    for i in range(segments.shape[0]):
        node_scores[i] = segments[i, 4] - (segments[i, 2] * ins_cost)

    pred.fill(-1)
    nb_node_scores.fill(-1e6)  # Must set to large negative, otherwise a value of zero can imply a path to that node

    best_score = 0  # Declare here in case only 1 alignment
    next_best_score = 0
    best_normal_orientation = 0  # Keep track of the best normal pairing score, F first R second
    normal_end_index = -1  # Keep track of the last normal-index for updating the normal-score later on

    # start from segment two because the first has been scored
    for i in range(1, segments.shape[0]):
        chr1 = segments[i, 0]
        pos1 = segments[i, 1]
        start1 = segments[i, 2]
        end1 = segments[i, 3]
        score1 = segments[i, 4]
        row_index1 = segments[i, 5]
        strand1 = segments[i, 6]
        r1 = segments[i, 7]

        p = -1  # -1 means there is no predecessor
        best_score = score1 - (start1 * ins_cost)  # Implies all preceding alignments skipped!
        next_best_score = - (start1 * ins_cost)  # Worst case

        # Walking backwards
        for j in range(i-1, -1, -1):

            chr2 = segments[j, 0]
            pos2 = segments[j, 1]
            start2 = segments[j, 2]
            end2 = segments[j, 3]
            score2 = segments[j, 4]
            row_index2 = segments[j, 5]
            strand2 = segments[j, 6]
            r2 = segments[j, 7]

            # Allow alignments with minimum sequence and max overlap
            if start1 > end2 - max_homology and end1 > end2 + min_aln and start1 - start2 > min_aln:

                if start1 > end2 and start1 - end2 > max_insertion:
                    continue

                # Microhomology and insertion lengths between alignments on same read only
                micro_h = 0
                ins = 0
                #if r1 == r2:
                micro_h = end2 - start1
                if micro_h < 0:
                    ins = fabs(micro_h)
                    micro_h = 0

                # Define jump cos
                proper_pair = 0
                if chr1 == chr2:
                    if r1 == r2:  # If the jump occurs on the same read, means there is an SV
                        jump_cost = U

                    else:
                        distance = fabs(pos1 - pos2)
                        proper_pair = is_proper_pair(distance, pos1, pos2, strand1, strand2, mu, sigma)
                        if proper_pair:
                            jump_cost = bwa_pair_score(distance, mu, sigma, match_score, U)

                            normal_jumps.add((i, j))  # Keep track of normal pairings on the path
                        else:
                            jump_cost = U

                else:
                    jump_cost = inter_cost + U

                # Calculate score, last_score = node_scores[j]
                total_cost = (micro_h * hom_cost) + (ins * ins_cost) + jump_cost
                S = score1 - total_cost  # Score 1 is 'ahead' in the sort order from sc2

                current_score = node_scores[j] + S

                # Update best score and next best score
                if current_score > best_score:
                    next_best_score = best_score
                    best_score = current_score

                    p = j

                elif current_score > next_best_score:
                    next_best_score = current_score

                if proper_pair:
                    if current_score > best_normal_orientation:
                        best_normal_orientation = current_score
                        normal_end_index = i  # index i is towards the right-hand side of the array

        node_scores[i] = best_score
        pred[i] = p
        nb_node_scores[i] = next_best_score

    # Update the score for jumping to the end of the sequence
    # Basically calculates what the best and secondary scores are for the 'virtual' end node
    cdef float node_to_end_cost
    cdef float secondary
    cdef float right_clip = 0
    cdef float cst = 0
    path_score = -(contig_length * ins_cost)  # Worst case
    secondary = float(path_score)
    end_i = -1
    for i in range(segments.shape[0]):

        node_to_end_cost = node_scores[i] - (ins_cost * (contig_length - segments[i, 3]))
        node_scores[i] = node_to_end_cost  # Need to know for secondary path
        if node_to_end_cost > path_score:
            path_score = node_to_end_cost
            end_i = i
        elif node_to_end_cost > secondary:
            secondary = node_to_end_cost

        if i == normal_end_index:  # Update the best normal score, deals with right-hand-side soft-clips
            best_normal_orientation = node_to_end_cost

    # Need to check if any branches from main path have a higher secondary than the 'virtual' end node
    # This can happen if the secondary path joins the main path within the graph, rather than at the end node.

    cdef float dis_2_end = 0  # Can be negative if the last node has an insertion before the end (soft-clip)
    cdef float s, potential_secondary

    # Get the path using traceback
    cdef vector[int] v
    cdef int normal_pairings = 0
    cdef int last_i
    cdef int next_i
    if end_i != -1:

        v.push_back(end_i)
        while True:
            last_i = v.back()
            next_i = pred[last_i]
            if next_i == -1:
                break
            v.push_back(next_i)

            if (last_i, next_i) in normal_jumps:
                normal_pairings += 1

            s = nb_node_scores[next_i]
            dis_2_end = path_score - node_scores[next_i]  # Distance of main path to the end
            potential_secondary = dis_2_end + s

            if potential_secondary > secondary:
                secondary = potential_secondary

    cdef np.ndarray[np.int_t, ndim=1] a = np.empty(len(v), dtype=np.int)
    for i in range(<int>v.size()):
        a[v.size() - 1 - i] = v[i]  # Virtual reversal of the indexes array

    if secondary < 0:
        secondary = 0
    if path_score < 0:
        path_score = 0
    if best_normal_orientation < 0:
        best_normal_orientation = 0
    # Return the index column
    return segments[a, 5], path_score, secondary, best_normal_orientation, normal_pairings


cpdef tuple process(dict rt):
    """
    Assumes that the reads are ordered read1 then read2, in the FR direction
    :param rt: Read_template object, contains all parameters within the pairing_params array
    :returns 5-tuple: list path, float length, float second_best, float dis_to_normal, float norm_pairings
    """

    pp = [float(i) for i in  rt["pairing_params"]]

    r1_len = rt['read1_length']
    r2_len = rt['read2_length']
    cdef int single_end = 0
    cdef int contig_l = 0

    cdef float max_insertion = pp[0]
    cdef float min_aln = pp[1]
    cdef float max_homology = pp[2]
    cdef float ins_cost = pp[3]
    cdef float hom_cost = pp[4]
    cdef float inter_cost = pp[5]
    cdef float U = pp[6]
    cdef float match_score = rt["match_score"]

    cdef np.ndarray[DTYPE_t, ndim=2] read1_arr, read2_arr, t, t2

    cdef float mu, sigma
    mu, sigma = rt['isize']

    t = rt['data'][:, 0:8]

    if not rt["paired_end"]:

        if rt["read1_unmapped"]:
            return np.array([0]).astype(np.float), 0, 0, 125, 0

        single_end = 1
        contig_l = r1_len
        return optimal_path(t, contig_l, mu, sigma, max_insertion, min_aln, max_homology, ins_cost,
                            hom_cost, inter_cost, U, match_score)

    else:
        if (r1_len is None and r2_len is None) or (rt["read1_unmapped"] and rt["read2_unmapped"]):
            return np.array([0, rt['first_read2_index']]).astype(np.float), 0, 0, 250, 0

        elif rt["read2_unmapped"]:

            single_end = 1
            contig_l = r1_len

            read1_arr = t[t[:, 7] == 1]

            path1, length1, second_best1, dis_to_normal1, norm_pairings1 = optimal_path(read1_arr, contig_l, mu, sigma,
                                                                                          max_insertion, min_aln,
                                                                                          max_homology, ins_cost,
                                                                                          hom_cost, inter_cost, U,
                                                                                          match_score)

            path1 = np.array(list(path1) + [rt['first_read2_index']]).astype(np.float)
            return path1, length1, second_best1, dis_to_normal1, norm_pairings1


        elif rt["read1_unmapped"]:

            single_end = 1
            contig_l = r2_len

            read2_arr = t[t[:, 7] == 2]

            path2, length2, second_best2, dis_to_normal2, norm_pairings2 = optimal_path(read2_arr, contig_l, mu, sigma,
                                                                                          max_insertion, min_aln,
                                                                                          max_homology, ins_cost,
                                                                                          hom_cost, inter_cost, U,
                                                                                          match_score)

            first_i = rt['first_read2_index']
            path2 = np.array([0] + list(path2)).astype(np.float)
            return path2, length2, second_best2, dis_to_normal2, norm_pairings2

        else:
            single_end = 0
            contig_l = r1_len + r2_len

            path1, length1, second_best1, dis_to_normal1, norm_pairings1 = optimal_path(t, contig_l, mu, sigma,
                                                                                              max_insertion, min_aln,
                                                                                              max_homology, ins_cost,
                                                                                              hom_cost, inter_cost, U,
                                                                                              match_score)


            if length1 == norm_pairings1 and len(path1) == 2 and length1 - second_best1 > U:
                return path1, length1, second_best1, dis_to_normal1, norm_pairings1


            read1_arr, read2_arr = [t[t[:, 7] == 1], t[t[:, 7] == 2]]
            read2_arr[:, 2:4] -= r1_len
            read1_arr[:, 2:4] += r2_len
            t2 = np.concatenate([read2_arr, read1_arr])

            path2, length2, second_best2, dis_to_normal2, norm_pairings2 = optimal_path(t2, contig_l, mu, sigma,
                                                                                              max_insertion, min_aln,
                                                                                              max_homology, ins_cost,
                                                                                              hom_cost, inter_cost, U,
                                                                                              match_score)

            if length2 > length1:
                return path2, length2, second_best2, dis_to_normal2, norm_pairings2
            else:
                return path1, length1, second_best1, dis_to_normal1, norm_pairings1
