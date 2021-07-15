import numpy as np
import pandas as pd
import click
import glob
from io import StringIO
from sys import stderr


@click.command()
@click.option("-w", help="Previous working directory for dysgu, converts all .bin files within -w", required=False, type=click.Path())
@click.option("-g", help="Converts all .bin files specified using pattern, e.g. 'wd/chr1*.bin", required=False, type=str)
@click.option("-b", help="Converts single .bin file", required=False, type=click.Path())
@click.option("--out-bin-size", help="Output bin size in base-pairs, must be a multiple of 100", type=click.IntRange(100, 10000000, clamp=True), default=100, show_default=True)
@click.option("--opp", help="When merging bins, apply this operation", type=click.Choice(["mean", "min", "max", "median"]), default="median", show_default=True)
@click.option("--sep", help="Separator", type=click.Choice([",", '\t', " "]), default="\t", show_default=True)
def convert2bed(w, g, b, out_bin_size, opp, sep):
    """Convert .dysgu_chrom.bin files to a bed file. Outputs '#chrom, start, end, coverage' columns to stdout"""
    if w is None and b is None and g is None:
        raise ValueError("Specify --wd or --bin")
    if out_bin_size % 100 > 0:
        raise ValueError("out-bin-size must be a multiple of 100 bp")
    if w is not None:
        bin_files = sorted(glob.glob(w + "/*dysgu_chrom.bin"))
    elif g is not None:
        bin_files = sorted(glob.glob(g))
    else:
        bin_files = [b]
    if len(bin_files) == 0:
        raise ValueError("No .bin files detected")

    block_size = int(out_bin_size / 10)
    print(f"Aggregating {block_size} consecutive bins", file=stderr)
    chroms = []
    starts = []
    ends = []
    vals = []
    for pth in bin_files:
        chrom_name = pth.split("/")[-1].split(".")[0]
        a = np.fromfile(pth, dtype="int16")
        for start in range(0, len(a), block_size):
            end = start + block_size
            if opp == "median":
                v = np.median(a[start:end])
            elif opp == "max":
                v = np.max(a[start:end])
            elif opp == "min":
                v = np.min(a[start:end])
            elif opp == "mean":
                v = np.mean(a[start:end])
            chroms.append(chrom_name)
            starts.append(start * 10)
            ends.append((start * 10) + out_bin_size)
            vals.append(v)

    df = pd.DataFrame()
    df["#chrom"] = chroms
    df["start"] = starts
    df["end"] = ends
    df["coverage"] = vals

    output = StringIO()
    df.to_csv(output, index=False, sep=sep)
    print(output.getvalue())


if __name__ == "__main__":
    convert2bed()
