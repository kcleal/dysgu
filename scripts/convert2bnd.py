from sys import stdout
from dysgu.view import read_from_inputfile
import click
import pysam
import copy


@click.command()
@click.option('-t', help="Comma separated list of SVTYPE's to convert to BND format", required=False,
              type=str, default='TRA', show_default=True)
@click.option('-o', help='Write vcf to output file', required=False, type=str, default='stdout', show_default=True)
@click.argument("reference", required=True, type=str)
@click.argument("input_vcf", required=False, type=str)
def conver2bnd(t, o, reference, input_vcf):
    if o == '-' or o == 'stdout':
        outf = stdout
    else:
        outf = open(o, 'w')

    fin = read_from_inputfile(input_vcf)
    targets = set(t.split(','))
    ref = pysam.FastaFile(reference)

    for line in fin:
        if line[0] == "#":
            outf.write(line)
        else:
            split = False
            for svt in targets:
                if svt in line:
                    split = svt
                    break
            if not split:
                outf.write(line)
                continue

            l = line.split('\t', 8)
            info = dict(tuple(i.split('=')) for i in l[7].split(';'))
            ct = info['CT']

            chrom2 = info['CHR2']
            pos2 = info["CHR2_POS"] if "CHR2_POS" in info else info["END"]
            p = f'{chrom2}:{pos2}'
            p2 = f'{l[0]}:{l[1]}'

            ref_base = l[3]

            ref_base2 = ref.fetch(chrom2, int(pos2), int(pos2)+1).upper()
            alt_bases = l[4]
            if alt_bases[0] != '<':
                add_ins = True
            else:
                add_ins = False
            if ct == '3to5':
                template = "{t}[{p}["
                template2 = "]{p}]{t}"
                if add_ins:
                    ref_base += alt_bases
                    ref_base2 = alt_bases + ref_base2

            elif ct == '5to3':
                template = "]{p}]{t}"
                template2 = "{t}[{p}["
                if add_ins:
                    ref_base = alt_bases + ref_base
                    ref_base2 += alt_bases
            elif ct == '3to3':
                template = "{t}]{p}]"
                template2 = "[{p}[{t}"
                if add_ins:
                    ref_base += alt_bases
                    ref_base2 = alt_bases + ref_base2
            else:
                template = "[{p}[{t}"
                template2 = "{t}]{p}]"
                if add_ins:
                    ref_base = alt_bases + ref_base
                    ref_base2 += alt_bases

            info['SVTYPE'] = 'BND'
            info2 = copy.deepcopy(info)

            info2['CHR2'] = l[0]
            if 'CHR2_POS' in info:
                info2['CHR2_POS'] = l[1]
                info2['END'] = str(int(pos2) + 1)
            else:
                info['END'] = l[1]

            l2 = copy.deepcopy(l)
            bnd1 = template.format(p=p, t=ref_base)
            bnd2 = template2.format(p=p2, t=ref_base2)

            l2[0] = chrom2
            l2[1] = pos2

            l[4] = bnd1
            l2[4] = bnd2

            l[7] = ';'.join([f'{i}={j}' for i, j in info.items()])
            l2[7] = ';'.join([f'{i}={j}' for i, j in info2.items()])

            l[2] += '_u'
            l2[2] += '_v'

            l = '\t'.join(l)
            l2 = '\t'.join(l2)

            outf.write(l)
            outf.write(l2)

    outf.close()


if __name__ == '__main__':
    conver2bnd()