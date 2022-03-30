import unittest
import pysam
import os
from dysgu import DysguSV, load_dysgu_vcf, merge_dysgu_df, merge_intervals
from tempfile import NamedTemporaryFile


class TestConstruct(unittest.TestCase):
    """ Test construction"""
    def test_api(self):

        test = os.path.abspath(os.path.dirname(__file__))

        bam = pysam.AlignmentFile(test + "/small.bam")
        genome = pysam.FastaFile(test + "/ref.fa")

        dysgu = DysguSV(genome, bam, remap=False, drop_gaps=False, min_support=3)
        f_iter = bam.fetch(until_eof=True)
        df = dysgu(f_iter)

        a = [("chr1", 1, 10, 0), ("chr1", 9, 11, 1), ("chr1", 20, 30, 2)]
        print(merge_intervals(a))
        print(merge_intervals(a, add_indexes=True))

        passed = df[df['filter'] == 'PASS']
        assert len(passed) == 1

        tempfile = NamedTemporaryFile("w")
        path = tempfile.name

        with open(path, "w") as out:
            dysgu.to_vcf(passed, out)

        df = load_dysgu_vcf(path, drop_na_columns=False)

        dysgu.set_option({"diploid": False})
        df2 = dysgu.apply_model(df)
        assert len(df2) == 1

        return 0

    def test_merge(self):
        bam = pysam.AlignmentFile("/Users/kezcleal/Desktop/HG002.bam")
        bam2 = pysam.AlignmentFile("/Users/kezcleal/Desktop/NA12878.bwa.bam")
        genome = pysam.FastaFile("/Users/kezcleal/Desktop/ucsc.hg19.fasta")

        dysgu = DysguSV(genome, bam, "hg002")
        r = ("chr1", 8217619, 8218348)
        r2 = ("chr1", 10_000_000, 15_000_000)
        df = dysgu(bam.fetch(*r))
        print(df)

        dysgu2 = DysguSV(genome, bam2, "n12878")
        df2 = dysgu2(bam2.fetch(*r))
        print(df2)

        df_merge = merge_dysgu_df(df, df2)
        print(df_merge)
        print("DOne")


def main():
    unittest.main()


if __name__ == "__main__":
    unittest.main()
