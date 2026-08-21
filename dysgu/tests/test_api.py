import unittest
import pysam
import os
import re
import pandas as pd
from io import StringIO
from dysgu import DysguSV, load_dysgu_vcf, merge_dysgu_df, merge_intervals
from dysgu.map_set_utils import EventResult
from dysgu.post_call import get_hp_format
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

    def test_gq_is_integer(self):
        """GQ is a reserved VCF FORMAT field and must be declared as Integer."""
        test = os.path.abspath(os.path.dirname(__file__))
        repo_root = os.path.dirname(os.path.dirname(test))
        bam = pysam.AlignmentFile(test + "/small.bam")
        genome = pysam.FastaFile(test + "/ref.fa")
        dysgu = DysguSV(genome, bam, remap=False, drop_gaps=False, min_support=3, metrics=False)

        df = load_dysgu_vcf(repo_root + "/test.dysgu1.9.0.vcf", drop_na_columns=False)
        # Duplicate the single row so we can test int/float/missing GQ coercion
        df = pd.concat([df, df.copy(), df.copy()], ignore_index=True)
        df.loc[df.index[0], "GQ"] = 30
        df.loc[df.index[1], "GQ"] = 30.0
        df.loc[df.index[2], "GQ"] = "."

        out = StringIO()
        dysgu.to_vcf(df, out)
        vcf_text = out.getvalue()

        # Header must declare GQ as Integer
        self.assertIn('##FORMAT=<ID=GQ,Number=1,Type=Integer', vcf_text)
        self.assertNotIn('##FORMAT=<ID=GQ,Number=1,Type=Float', vcf_text)

        # All non-missing GQ values must be integer strings
        fmt_idx = None
        for line in vcf_text.splitlines():
            if line.startswith("#CHROM"):
                fmt_idx = line.split("\t").index("FORMAT")
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            fmt = fields[fmt_idx].split(":")
            sample = fields[fmt_idx + 1].split(":")
            gq = sample[fmt.index("GQ")]
            if gq != ".":
                self.assertTrue(re.fullmatch(r"-?\d+", gq),
                                f"GQ value is not an integer: {gq}")

        # Round-trip through load_dysgu_vcf keeps integer GQ
        with NamedTemporaryFile("w", suffix=".vcf", delete=False) as tmp:
            tmp.write(vcf_text)
            tmp_path = tmp.name
        df2 = load_dysgu_vcf(tmp_path, drop_na_columns=False)
        os.unlink(tmp_path)
        self.assertTrue(all(str(v) == "." or re.fullmatch(r"-?\d+", str(v))
                            for v in df2["GQ"]),
                        "GQ values became non-integer after round-trip")

    def test_hp_unphased_only(self):
        """When only unphased reads support an SV, HP must be '_N', never blank."""
        e = EventResult()
        e.GT = "0/1"
        e.haplotype_counts = {'u': 6}
        e.phase_set_counts = {}
        get_hp_format([e])
        self.assertEqual(e.haplotype, "_6",
                         "HP should be '_N' when only unphased reads are present")

        # Mixed phased + unphased
        e2 = EventResult()
        e2.GT = "0/1"
        e2.haplotype_counts = {1: 4, 'u': 2}
        e2.phase_set_counts = {}
        get_hp_format([e, e2])
        self.assertEqual(e2.haplotype, "4_2")


def main():
    unittest.main()


if __name__ == "__main__":
    unittest.main()
