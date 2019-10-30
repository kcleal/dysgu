#cython: language_level=3

"""
Find the best path through a collection of alignments
Works with paried-end reads rather or single contigs. Borrows the pairing heuristic from bwa.
"""
from __future__ import absolute_import
import sys
import numpy as np
import os
import click
from dysgu import align_path


def read_orientations(t, r1_len, r2_len):
    """
    # This code prooves that the orientation of the read dowesnt matter, only the order of the reads:
    # Permute all orientations, [AB|CD] + AB|DC, BA|DC, BA|CD
    ab, cd = [t[t[:, 7] == 1], t[t[:, 7] == 2]]

    # Make AB reversed
    ba = ab.copy()
    s, e = ba[:, 2].copy(), ba[:, 3].copy()
    ba[:, 2] = r1_l - e
    ba[:, 3] = r1_l - s
    ba = np.flipud(ba)

    # Make CD reversed
    dc = cd.copy()
    dc[:, 2:4] -= r1_l
    s, e = dc[:, 2].copy(), dc[:, 3].copy()
    dc[:, 2] = r2_l - e
    dc[:, 3] = r2_l - s
    dc[:, 2:4] += r1_l
    dc = np.flipud(dc)

    # Add combos
    ori.append(np.concatenate([ab, dc]))
    ori.append(np.concatenate([ba, dc]))
    ori.append(np.concatenate([ba, cd]))

    These 3 plus the original, produce two sets of identical paths; therefore only the order of reads is important
    """
    yield t

    read1_arr, read2_arr = [t[t[:, 7] == 1], t[t[:, 7] == 2]]
    read2_arr[:, 2:4] -= r1_len
    read1_arr[:, 2:4] += r2_len

    yield np.concatenate([read2_arr, read1_arr])


def process(rt):
    """
    Assumes that the reads are ordered read1 then read2, in the FR direction
    :param rt: Read_template object, contains all parameters within the pairing_params array
    """
    # if rt["name"] == "HISEQ2500-10:539:CAV68ANXX:7:2204:21140:9933":
    #     click.echo(rt, err=True)
    #     quit()

    r1_len = rt['read1_length']
    r2_len = rt['read2_length']
    if not rt["paired_end"]:
        single_end = True
        contig_l = r1_len
    else:
        if r2_len is None and r1_len:
            single_end = True
            contig_l = r1_len
        elif r1_len is None and r2_len:
            single_end = True
            contig_l = r2_len
        elif r1_len is None and r2_len is None:
            return False
        else:
            single_end = False
            contig_l = r1_len + r2_len

    mu, sigma = rt['isize']

    pp = [float(i) for i in  rt["pairing_params"]]
    max_insertion = pp[0]
    min_aln = pp[1]
    max_homology = pp[2]
    ins_cost = pp[3]
    hom_cost = pp[4]  # 2
    inter_cost = pp[5]  # 10
    U = pp[6]
    match_score = rt["match_score"]

    args = [contig_l, mu, sigma, max_insertion, min_aln, max_homology, ins_cost,
            hom_cost, inter_cost, U, match_score]

    table = rt['data'][:, range(8)]

    if not single_end:

        # If it is unclear which read comes first this function can be used to generate both orientations:
        both_ways = []
        for r in read_orientations(table, r1_len, r2_len):
            a_res = align_path.optimal_path(r, *args)
            if len(a_res) > 0:
                if a_res[1] == a_res[4] and len(a_res[0]) == 2 and a_res[1] - a_res[2] > U:
                    # Cant do better. Normal pairing with 2 good alignments
                    return a_res
            both_ways.append(a_res)

        if len(both_ways) == 0:
            return False
        path, length, second_best, dis_to_normal, norm_pairings = sorted(both_ways, key=lambda x: x[1])[-1]  # Best

    else:
        path, length, second_best, dis_to_normal, norm_pairings = align_path.optimal_path(table, *args)

    if int(length) < int(second_best):
        sys.stderr.write("WARNING: primary path < secondary path\n")

    # if length < 0:  # Todo make a minimum allowed alignment score
    #     return False  # Drop template, couldnt align properly

    # Todo second best can be negative?
    return path, length, second_best, dis_to_normal, norm_pairings


if __name__ == "__main__":
    array = np.array
    rt = {'paired_end': 1, 'bias': 1.15, 'fq_read2_seq': 0, 'isize': (250.0, 50.0), 'read2_q': 'BBCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGCGGGEGGFGGGGGGGGGGCFGE', 'max_d': 450.0, 'read2_seq': 'TCTGTAGATCTTCTTTGGTGAAGTGTCTGTTCAGATCTGTGTGCATTTTTAATTGACAAATGACACAGAAGGTTGATATACACAGAAGAAATGACAATGTGGATTTCTTAATATTTACAGTTTAT', 'chrom_ids': {'chr6': 10, 'chr4': 8, 'chr2': 4, 'chr1': 1, 'chr8': 9, 'chr13': 6, 'chr11': 13, 'chr10': 7, 'chr17': 3, 'chr16': 12, 'chr20': 11, 'chr21': 0, 'chr22': 2, 'chr19': 5}, 'read2_length': 125, 'passed': 0, 'replace_hard': 0, 'read2_reverse': 0, 'inputdata': [['65', 'chr21', '46696680', '0', '125M', 'chr19', '248503', '0', 'GATCTGTATTTTGCAAATATTTTCTTCAATATGTGGCTTGTCTTTTTGTTCTCTTGACAAAGTCTCTTCCAGAGTATAAACTGTAAATATTAAGAAATCCACATTGTCATTTCTTCTGTGTATAT', 'CBCCCGGGGGGGGGGGGGGGGGGGG@GGGGGGGGGGGGGGGGGGGGGGGGGGGGGFEGGGGGGGGGGG1<@FDGGEGGGG>FGGCG>GGGGFF0B>DG>>EDFGGGGD@GGGGGGEFGGGG>GDG', 'NM:i:1', 'MD:Z:60G64', 'AS:i:120', 'XS:i:120'], ['321', 'chr1', '248942783', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:1', 'MD:Z:60G64', 'AS:i:120'], ['321', 'chr22', '50805037', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:1', 'MD:Z:60G64', 'AS:i:120'], ['321', 'chr17', '83246855', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:2', 'MD:Z:41G18G64', 'AS:i:115'], ['321', 'chr2', '242180718', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:2', 'MD:Z:55A4G64', 'AS:i:115'], ['321', 'chr19', '58605076', '0', '125M', '=', '248503', '-58356574', '*', '*', 'NM:i:2', 'MD:Z:55A4G64', 'AS:i:115'], ['321', 'chr13', '114351467', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:2', 'MD:Z:55A4G64', 'AS:i:115'], ['321', 'chr10', '133784627', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:2', 'MD:Z:55A4G64', 'AS:i:115'], ['321', 'chr4', '190201925', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:2', 'MD:Z:55A4G64', 'AS:i:115'], ['337', 'chr19', '248523', '0', '125M', '=', '248503', '-145', '*', '*', 'NM:i:3', 'MD:Z:64C4T8C46', 'AS:i:110'], ['337', 'chr8', '208316', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:3', 'MD:Z:64C4T8C46', 'AS:i:110'], ['321', 'chr6', '170607571', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:3', 'MD:Z:55A4G21T42', 'AS:i:110'], ['337', 'chr2', '113605836', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:3', 'MD:Z:64C4T25G29', 'AS:i:110'], ['321', 'chr20', '64284388', '0', '125M', 'chr19', '248503', '0', '*', '*', 'NM:i:5', 'MD:Z:55A1G2G21T4G37', 'AS:i:100'], ['129', 'chr19', '248503', '0', '55S70M', 'chr21', '46696680', '0', 'TCTGTAGATCTTCTTTGGTGAAGTGTCTGTTCAGATCTGTGTGCATTTTTAATTGACAAATGACACAGAAGGTTGATATACACAGAAGAAATGACAATGTGGATTTCTTAATATTTACAGTTTAT', 'BBCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGCGGGEGGFGGGGGGGGGGCFGE', 'NM:i:0', 'MD:Z:70', 'AS:i:70', 'XS:i:70', 'SA:Z:chr21,46696564,+,55M70S,0,0;'], ['401', 'chr1', '248942858', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr22', '50805112', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr2', '242180793', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr21', '46696755', '0', '70M55S', '=', '46696680', '-145', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr4', '190202000', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr13', '114351542', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['385', 'chr2', '113605816', '0', '55S70M', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr17', '83246930', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['385', 'chr8', '208296', '0', '55S70M', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr19', '58605151', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr10', '133784702', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:70', 'AS:i:70'], ['401', 'chr6', '170607646', '0', '70M55S', 'chr21', '46696680', '0', '*', '*', 'NM:i:2', 'MD:Z:7T49T12', 'AS:i:60'], ['2177', 'chr21', '46696564', '0', '55M70S', '=', '46696680', '117', 'TCTGTAGATCTTCTTTGGTGAAGTGTCTGTTCAGATCTGTGTGCATTTTTAATTGACAAATGACACAGAAGGTTGATATACACAGAAGAAATGACAATGTGGATTTCTTAATATTTACAGTTTAT', 'BBCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGCGGGEGGFGGGGGGGGGGCFGE', 'NM:i:0', 'MD:Z:55', 'AS:i:55', 'XS:i:55', 'SA:Z:chr19,248503,+,55S70M,0,0;'], ['385', 'chr22', '50804924', '0', '55M70S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:55', 'AS:i:55'], ['385', 'chr1', '248942667', '0', '55M70S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:55', 'AS:i:55'], ['401', 'chr2', '113606022', '0', '70S55M', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:55', 'AS:i:55'], ['385', 'chr17', '83246739', '0', '55M70S', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:55', 'AS:i:55'], ['385', 'chr19', '58604960', '0', '55M70S', 'chr21', '46696680', '0', '*', '*', 'NM:i:1', 'MD:Z:21G33', 'AS:i:50'], ['385', 'chr20', '64284271', '0', '55M70S', 'chr21', '46696680', '0', '*', '*', 'NM:i:1', 'MD:Z:21G33', 'AS:i:50'], ['401', 'chr1', '77234119', '0', '91S34M', 'chr21', '46696680', '0', '*', '*', 'NM:i:0', 'MD:Z:34', 'AS:i:34'], ['401', 'chr16', '1251048', '0', '87S30M2I6M', 'chr21', '46696680', '0', '*', '*', 'NM:i:2', 'MD:Z:36', 'AS:i:31'], ['401', 'chr11', '126671620', '0', '87S38M', 'chr21', '46696680', '0', '*', '*', 'NM:i:2', 'MD:Z:31A0C5', 'AS:i:31'], ['385', 'chr16', '1235361', '0', '38M87S', 'chr21', '46696680', '0', '*', '*', 'NM:i:2', 'MD:Z:3C2T31', 'AS:i:31']], 'fq_read1_q': 0, 'fq_read2_q': 0, 'read1_reverse': 0, 'read1_q': 'CBCCCGGGGGGGGGGGGGGGGGGGG@GGGGGGGGGGGGGGGGGGGGGGGGGGGGGFEGGGGGGGGGGG1<@FDGGEGGGG>FGGCG>GGGGFF0B>DG>>EDFGGGGD@GGGGGGEFGGGG>GDG', 'read1_length': 125, 'data': array([[ 0.00000000e+00,  4.66966800e+07,  0.00000000e+00,
         1.25000000e+02,  1.37999997e+02,  0.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.20000000e+02],
       [ 1.00000000e+00,  2.48942783e+08,  0.00000000e+00,
         1.25000000e+02,  1.20000000e+02,  1.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.20000000e+02],
       [ 2.00000000e+00,  5.08050370e+07,  0.00000000e+00,
         1.25000000e+02,  1.20000000e+02,  2.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.20000000e+02],
       [ 3.00000000e+00,  8.32468550e+07,  0.00000000e+00,
         1.25000000e+02,  1.15000000e+02,  3.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.15000000e+02],
       [ 4.00000000e+00,  2.42180718e+08,  0.00000000e+00,
         1.25000000e+02,  1.15000000e+02,  4.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.15000000e+02],
       [ 5.00000000e+00,  5.86050760e+07,  0.00000000e+00,
         1.25000000e+02,  1.15000000e+02,  5.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.15000000e+02],
       [ 6.00000000e+00,  1.14351467e+08,  0.00000000e+00,
         1.25000000e+02,  1.15000000e+02,  6.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.15000000e+02],
       [ 7.00000000e+00,  1.33784627e+08,  0.00000000e+00,
         1.25000000e+02,  1.15000000e+02,  7.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.15000000e+02],
       [ 8.00000000e+00,  1.90201925e+08,  0.00000000e+00,
         1.25000000e+02,  1.15000000e+02,  8.00000000e+00,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.15000000e+02],
       [ 5.00000000e+00,  2.48523000e+05,  0.00000000e+00,
         1.25000000e+02,  1.10000000e+02,  9.00000000e+00,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.10000000e+02],
       [ 9.00000000e+00,  2.08316000e+05,  0.00000000e+00,
         1.25000000e+02,  1.10000000e+02,  1.00000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.10000000e+02],
       [ 1.00000000e+01,  1.70607571e+08,  0.00000000e+00,
         1.25000000e+02,  1.10000000e+02,  1.10000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.10000000e+02],
       [ 4.00000000e+00,  1.13605836e+08,  0.00000000e+00,
         1.25000000e+02,  1.10000000e+02,  1.20000000e+01,
        -1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.10000000e+02],
       [ 1.10000000e+01,  6.42843880e+07,  0.00000000e+00,
         1.25000000e+02,  1.00000000e+02,  1.30000000e+01,
         1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
         1.00000000e+02],
       [ 0.00000000e+00,  4.66965640e+07,  1.25000000e+02,
         1.80000000e+02,  6.32499987e+01,  2.70000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.50000000e+01],
       [ 2.00000000e+00,  5.08049240e+07,  1.25000000e+02,
         1.80000000e+02,  5.50000000e+01,  2.80000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.50000000e+01],
       [ 1.00000000e+00,  2.48942667e+08,  1.25000000e+02,
         1.80000000e+02,  5.50000000e+01,  2.90000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.50000000e+01],
       [ 4.00000000e+00,  1.13606022e+08,  1.25000000e+02,
         1.80000000e+02,  5.50000000e+01,  3.00000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.50000000e+01],
       [ 3.00000000e+00,  8.32467390e+07,  1.25000000e+02,
         1.80000000e+02,  5.50000000e+01,  3.10000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.50000000e+01],
       [ 5.00000000e+00,  5.86049600e+07,  1.25000000e+02,
         1.80000000e+02,  5.00000000e+01,  3.20000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.00000000e+01],
       [ 1.10000000e+01,  6.42842710e+07,  1.25000000e+02,
         1.80000000e+02,  5.00000000e+01,  3.30000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         5.00000000e+01],
       [ 1.00000000e+00,  7.72341190e+07,  1.25000000e+02,
         1.59000000e+02,  3.40000000e+01,  3.40000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         3.40000000e+01],
       [ 1.20000000e+01,  1.25104800e+06,  1.25000000e+02,
         1.63000000e+02,  3.10000000e+01,  3.50000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         3.10000000e+01],
       [ 1.30000000e+01,  1.26671620e+08,  1.25000000e+02,
         1.63000000e+02,  3.10000000e+01,  3.60000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         3.10000000e+01],
       [ 1.20000000e+01,  1.23536100e+06,  1.25000000e+02,
         1.63000000e+02,  3.10000000e+01,  3.70000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         3.10000000e+01],
       [ 0.00000000e+00,  4.66967550e+07,  1.80000000e+02,
         2.50000000e+02,  8.04999983e+01,  1.80000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 5.00000000e+00,  2.48503000e+05,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  1.40000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 1.00000000e+00,  2.48942858e+08,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  1.50000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 2.00000000e+00,  5.08051120e+07,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  1.60000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 4.00000000e+00,  2.42180793e+08,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  1.70000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 8.00000000e+00,  1.90202000e+08,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  1.90000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 6.00000000e+00,  1.14351542e+08,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  2.00000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 4.00000000e+00,  1.13605816e+08,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  2.10000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 3.00000000e+00,  8.32469300e+07,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  2.20000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 9.00000000e+00,  2.08296000e+05,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  2.30000000e+01,
         1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 5.00000000e+00,  5.86051510e+07,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  2.40000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 7.00000000e+00,  1.33784702e+08,  1.80000000e+02,
         2.50000000e+02,  7.00000000e+01,  2.50000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         7.00000000e+01],
       [ 1.00000000e+01,  1.70607646e+08,  1.80000000e+02,
         2.50000000e+02,  6.00000000e+01,  2.60000000e+01,
        -1.00000000e+00,  2.00000000e+00,  0.00000000e+00,
         6.00000000e+01]]), 'name': 'HISEQ2500-10:539:CAV68ANXX:7:2204:21140:9933', 'fq_read1_seq': 0, 'match_score': 1.0, 'read1_seq': 'GATCTGTATTTTGCAAATATTTTCTTCAATATGTGGCTTGTCTTTTTGTTCTCTTGACAAAGTCTCTTCCAGAGTATAAACTGTAAATATTAAGAAATCCACATTGTCATTTCTTCTGTGTATAT', 'last_seen_chrom': 'chr16', 'score_mat': {}, 'pairing_params': (150.0, 17.0, 150.0, 0.1, 3.0, 2.0, 9.0)}

    print(process(rt))

    for row in rt["data"]:
        print(list(row.astype(int)))

