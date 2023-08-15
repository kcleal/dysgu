import click
from sys import stderr
import pysam
from dysgu.coverage import index_stats


@click.command()
@click.argument('alignment_file', required=True, type=click.Path(exists=False))
@click.option("-y", "--y", help="Max-cov is estimated as {mean-coverage} * y", required=False, type=float,
              show_default=True, default=6)
def suggest_max_coverage(alignment_file, y):
    """Estimate a max-coverage value for use with dysgu. Mean genome coverage is estimated from the index file, so
    will only be useful for whole-genome alignment files"""
    f = pysam.AlignmentFile(alignment_file)
    cov, read_length = index_stats(f)
    cov = round(cov, 2)
    read_length = round(read_length, 1)
    max_cov = round(cov * y)
    print(f"Read-length {read_length} bp, mean whole-genome coverage estimate: {cov}, max-cov ~ {max_cov}", file=stderr)
    print(max_cov)
    return max_cov


if __name__ == "__main__":
    suggest_max_coverage()
