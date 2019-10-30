from __future__ import absolute_import
import click
from click.testing import CliRunner
import datetime
import os
import time
from multiprocessing import cpu_count
from subprocess import Popen, PIPE, check_call
from dysgu import find_pairs
from dysgu import input_stream_alignments, data_io, cluster, get_reads, view
import pkg_resources
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

cpu_range = click.IntRange(min=1, max=cpu_count())

defaults = {
            "clip_length": 21,
            "mapper": "bwamem",
            "map_script": None,
            "procs": 1,
            "dest": None,
            "post_fix": "dysgu",
            "search": None,
            "exclude": None,
            "include": None,
            "paired": "True",
            "read_length": 125.,
            "max_insertion": 150.,
            "min_aln": 17.,
            "max_overlap": 150.,
            "ins_cost": 0.1,
            "ol_cost": 1.,
            "inter_cost": 2.,
            "u": 9.,
            "match_score": 1.,
            "bias": 1.15,
            "output": "-",
            "svs_out": "-",
            "replace_hardclips": "False",
            "fq1": None,
            "fq2": None,
            "max_cov": 150,
            "buffer_size": 500000,
            "min_support": 3,
            "I": "210,175",
            "mark_dups": "True",
            "model": None,
            "mq": None,
            "max_tlen": 800,
            "template_size": "",
            "z_depth": 2,
            "z_breadth": 3,
            "regions_only": "False",
            "soft_search": "True"
            }

align_args = {}

version = pkg_resources.require("dysgu")[0].version

# Todo Make vcf output option


def pipeline(kwargs):
    t0 = time.time()
    click.echo("Running dysgu pipeline", err=True)
    if kwargs["bam"] is None:
        raise IOError("Error: Input bam is None")

    if not os.path.exists(kwargs["bam"] + ".bai"):
        raise IOError("Input .bai index file not found.")

    data_io.mk_dest(kwargs["dest"])

    #other_kwargs = find_pairs.process(kwargs)
    other_kwargs = get_reads.process(kwargs)
    kwargs.update(other_kwargs)

    single = True if kwargs["procs"] == 1 else False

    process, used_procs, remaining_procs = launch_external_mapper(kwargs)
    kwargs["procs"] = remaining_procs

    if kwargs["mapper"] == "bwamem" and kwargs["map_script"] is None:
        kwargs["fq1"] = None
        kwargs["fq2"] = None
    else:
        kwargs["fq1"] = kwargs["out_pfix"] + "1.fq"
        kwargs["fq2"] = kwargs["out_pfix"] + "2.fq"

    kwargs["sam"] = iter(process.stdout.readline, "")
    kwargs["output"] = kwargs["out_pfix"] + ".sam"

    input_stream_alignments.process_reads(kwargs)

    process.kill()

    if not single:
        kwargs["procs"] = remaining_procs + used_procs - 1

    sort_and_index(kwargs)

    # Clean up
    if kwargs["fq1"] is not None:
        os.remove(kwargs["fq1"])
        os.remove(kwargs["fq2"])

    if kwargs["mapper"] == "last":
        os.remove(kwargs["out_pfix"] + ".dict")

    kwargs["sv_aligns"] = kwargs["out_pfix"] + ".srt.bam"
    kwargs["raw_aligns"] = kwargs["bam"]

    # cluster_old.cluster_reads(kwargs)
    click.echo("dysgu run {} completed in {} h:m:s\n".format(kwargs["bam"],
                                                            str(datetime.timedelta(seconds=int(time.time() - t0)))),
               err=True)


def launch_external_mapper(kwargs):
    """Run a shell script to map interleaved .fastq files to .sam format. Uses a basic bwa-mem by default.
    A custom shell script can be provided but must take positional arguments as:
    $1 reference genome
    $2 .fastq file; interleaved if paired-end, otherwise single end reads"""

    # if kwargs["procs"] > 1:
    #     p = int(kwargs["procs"]) - 1
    #     other_p = 1  # kwargs["procs"] - p
    # else:
    #     p = kwargs["procs"]
    #     other_p = 1

    p = int(kwargs["procs"])
    if not kwargs["map_script"] and kwargs["mapper"] == "bwamem":
        command = "bwa mem -Y -t {procs} -a {ref} {s}1.fq {s}2.fq".format(procs=p - 1 if p > 1 else 1,
                                                                          ref=kwargs["reference"],
                                                                          s=kwargs["fastq"])

    elif not kwargs["map_script"] and kwargs["mapper"] == "last":

        command = "fastq-interleave {s}1.fq {s}2.fq \
        | lastal -k2 -l11 -Q1 -D10000 -K8 -C8 -i10M -r1 -q4 -a6 -b1 -P{procs} {ref} \
        | last-map-probs -m 1 -s 1 | maf-convert -f {d}.dict sam".format(procs=p,
                                                                         ref=kwargs["reference"],
                                                                         d=kwargs["out_pfix"],
                                                                         s=kwargs["fastq"])

    else:
        command = "bash {script} {ref} {s}1.fq {s}2.fq {procs}".format(script=kwargs["map_script"],
                                                                       procs=p,
                                                                       ref=kwargs["reference"],
                                                                       s=kwargs["fastq"])

    click.echo("Mapping command:\n" + command, err=True)

    proc = Popen(command, stdout=PIPE, shell=True, bufsize=0, universal_newlines=True)

    return proc, p, 1  # Note dysgu align always has 1 core assigned


def sort_and_index(kwargs):
    """Convenience function to sort and index a sam file, then remove the input sam file"""
    c = "samtools view -Sh {fix}.sam | \
    samtools sort -@ {p} -o {fix}.srt.bam - ; \
    samtools index -@ {p} {fix}.srt.bam"
    c = c.format(fix=kwargs["out_pfix"], p=kwargs["procs"])
    click.echo(c, err=True)
    check_call(c, shell=True)
    os.remove(kwargs["output"])


# User Interface:
# ----------------------------------------------------------------------------------------------------------------------
def apply_ctx(ctx, kwargs):
    click.echo("[dysgu] Version: {}".format(version), err=True)
    ctx.ensure_object(dict)
    if len(ctx.obj) == 0:  # When run is invoked from cmd line, else run was invoked from test function
        for k, v in list(defaults.items()) + list(kwargs.items()):
            ctx.obj[k] = v
    i, j = map(float, ctx.obj["I"].split(","))
    if "insert_median" not in ctx.obj:
        ctx.obj["insert_median"] = i
        ctx.obj["insert_stdev"] = j

    return ctx


@click.group(chain=False, invoke_without_command=False)
@click.version_option()
def cli():
    """Dysgu-SV is a set of tools for mapping and calling structural variants from sam/bam/cram files"""
    pass


# @cli.command("run")
# @click.argument('reference', required=True, type=click.Path(exists=False))
# @click.argument("bam", required=True, type=click.Path(exists=True))
# @click.option("-o", "--svs-out", help="Structural variants output, default=stdout", required=False, type=click.Path())
# @click.option('--include', help=".bed file, limit calls to regions", default=None, type=click.Path(exists=True))
# @click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
# @click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
#               default=None, type=click.Path(exists=True))
# @click.option('--clip-length', help="Minimum soft-clip length; >= threshold are kept", default=defaults["clip_length"],
#               type=int, show_default=True)
# @click.option('--mapper', help="External mapper to use for re-alignment", default=defaults["mapper"],
#               type=click.Choice(['bwamem', 'last']), show_default=True)
# @click.option('--map-script', help="""External shell script for mappping fastq files. \
# Overrides --mapper argument. Script must take positional arguments as: $1 reference genome; \
# $2 read1.fastq file $3 read2.fastq if paired-end \
# $4 threads to use""", default=None, type=click.Path(exists=True))
# @click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
# @click.option('--dest', help="Destination folder to use/create for saving results. Defaults to current directory",
#               default=None, type=click.Path())
# @click.option('-I', help="Insert size and stdev as 'FLOAT,FLOAT'. If not provided, automatically inferred",
#               default=defaults["I"], type=str)
# @click.option("--model", help="A model trained with dysgu train", default=defaults["model"],
#               type=click.Path(), show_default=True)
# @click.option('--mq', help="MapQ recalibration model", default=defaults["mq"],
#               type=click.Path(), show_default=True, required=False)
# @click.pass_context
# def run_command(ctx, **kwargs):
#     """Run the dysgu pipeline."""
#     ctx = apply_ctx(ctx, kwargs)
#     pipeline(ctx.obj)


@cli.command("reads")
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.option('--post-fix', help="Post-fix to use", default='dysgu', type=str, show_default=True)
@click.option('--clip-length', help="Minimum soft-clip length, > threshold are kept", default=defaults["clip_length"],
              type=int, show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=defaults["procs"], show_default=True)
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
              default=None, type=click.Path(exists=True))
@click.option("--soft-search",  help="Set to True will collect all reads with one or more alignments in\
 --search regions",
              default=defaults["soft_search"], type=click.Choice(["True", "False"]), show_default=True)
# @click.option('--bam-only', help="Output bam instead of fastq", default="False", type=click.Choice(["True", "False"]),
#               show_default=True)
@click.option('--dest', help="Destination folder to use/create. Defaults to current directory",
              default=None, type=click.Path())
@click.pass_context
def find_reads(ctx, **kwargs):
    """Finds read templates that are discordant or have an alignment with a soft-clip > '--clip-length'"""
    # Add arguments to context insert_median, insert_stdev, read_length, out_name
    ctx = apply_ctx(ctx, kwargs)
    return find_pairs.process(ctx.obj)


@cli.command("sv2fq")
@click.argument('bam', required=True, type=click.Path(exists=True))
@click.option('--post-fix', help="Post fix to tag temp files with", default='dysgu', type=str, show_default=True)
@click.option('--clip-length', help="Minimum soft-clip length, > threshold are kept", default=defaults["clip_length"],
              type=int, show_default=True)
@click.option("--paired", help="Paired end reads (or single)", default=defaults["paired"],
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("-f", "--out-format", help="Output format", default="fq",
              type=click.Choice(["fq", "fasta"]), show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=defaults["procs"], show_default=True)
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
              default=None, type=click.Path(exists=True))
@click.option('--dest', help="Destination folder to use/create for saving results. Defaults to current directory",
              default=None, type=click.Path())
@click.pass_context
def sv2fq(ctx, **kwargs):
    """Filters input .bam/.cram for read-pairs that are discordant or have a soft-clip of length > '--clip-length',
    writes two .fq files, unless --paired False"""
    # Add arguments to context insert_median, insert_stdev, read_length, out_name
    ctx = apply_ctx(ctx, kwargs)
    return get_reads.process(ctx.obj)


@cli.command("choose")
@click.argument("sam", type=click.File('r', encoding="ascii"), required=True)
@click.argument("output", required=False, type=click.Path())
@click.option("--paired", help="Paired end reads (or single)", default=defaults["paired"],
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-I', help="Insert size and stdev as 'FLOAT,FLOAT'",
              default=defaults["I"], type=str, show_default=True)
@click.option('--mq', help="MapQ recalibration model", default=defaults["mq"],
              type=click.Path(), show_default=True, required=False)
@click.option("--replace-hardclips",  help="Replace hard-clips with soft-clips when possible",
              default=defaults["replace_hardclips"], type=click.Choice(["True", "False"]), show_default=True)
@click.option("--fq1",  help="Fastq/fasta reads 1, used to add soft-clips to all hard-clipped read 1 alignments",
              default=defaults["fq1"], type=click.Path(), show_default=True)
@click.option("--fq2",  help="Fastq/fasta reads 2, used to add soft-clips to all hard-clipped read 2 alignments",
              default=defaults["fq2"], type=click.Path(), show_default=True)
@click.option("--max-insertion", help="Maximum insertion within read", default=defaults["max_insertion"], type=float,
              show_default=True)
@click.option("--min-aln", help="Minimum alignment length", default=defaults["min_aln"], type=float, show_default=True)
@click.option("--max-overlap", help="Maximum overlap between successive alignments", default=defaults["max_overlap"],
              type=float, show_default=True)
@click.option("--ins-cost", help="Insertion cost", default=defaults["ins_cost"], type=float, show_default=True)
@click.option("--ol-cost", help="Overlapping alignment cost", default=defaults["ol_cost"], type=float)
@click.option("--inter-cost", help="Cost of inter-chromosomal jump", default=defaults["inter_cost"], type=float,
              show_default=True)
@click.option("--u", help="Pairing cost", default=defaults["u"], type=float, show_default=True)
@click.option("--match-score", help="Matched base score", default=defaults["match_score"],
              type=float, show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option('--include', help=".bed file, elevate alignment scores in these regions. Determined by '--bias'",
              default=None, type=click.Path(exists=True))
@click.option("--bias", help="""Multiply match score by bias if alignment falls within regions .bed file
Unused if .bed not provided""", default=defaults["bias"], type=float, show_default=True)
@click.pass_context
def dysgu_aligner(ctx, **kwargs):
    """Choose an optimal set of alignments from from a collection of candidate alignments.
    If reads are paired, alignments must be sorted by read-name with the bit flag
    designating read_1 vs read_2."""
    ctx = apply_ctx(ctx, kwargs)
    input_stream_alignments.process_reads(ctx.obj)


@cli.command("call")
@click.argument('sv-aligns', required=True, type=click.Path(exists=True))
@click.option("-o", "--svs-out", help="Output file, [default: stdout]", required=False, type=click.Path())
@click.option('--clip-length', help="Minimum soft-clip length, > threshold are kept.", default=defaults["clip_length"],
              type=int, show_default=True)
@click.option('--max-cov', help="Regions with > max-cov that do no overlap 'include' are discarded",
              default=defaults["max_cov"], type=float, show_default=True)
@click.option('--max-tlen', help="Maximum template length to consider when calculating paired-end template size",
              default=defaults["max_tlen"], type=int, show_default=True)
@click.option('--min-support', help="Minimum number of reads per SV",
              default=defaults["min_support"], type=int, show_default=True)
@click.option('--z-depth', help="Minimum minimizer depth across alignments",
              default=defaults["z_depth"], type=int, show_default=True)
@click.option('--z-breadth', help="Minimum number of minimizers shared between a pair of alignments",
              default=defaults["z_breadth"], type=int, show_default=True)
@click.option("--template-size", help="Manually set insert size and stdev as 'INT,INT'",
              default=defaults["template_size"], type=str, show_default=False)
@click.option('--regions-only', help="If --include is provided, call only events within target regions",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=defaults["procs"], show_default=True)
@click.option('--include', help=".bed file, limit calls to regions", default=None, type=click.Path(exists=True))
@click.option('--dest', help="Folder to use/create for saving results. Defaults to current directory",
              default=None, type=click.Path())
@click.option("--buffer-size", help="Number of alignments to buffer", default=defaults["buffer_size"],
              type=int, show_default=True)
@click.option("--model", help="A model trained with dysgu train", default=defaults["model"],
              type=click.Path(), show_default=True)
@click.pass_context
def call_events(ctx, **kwargs):
    """Call structural vaiants"""
    # Create dest in not done so already
    ctx = apply_ctx(ctx, kwargs)
    cluster.cluster_reads(ctx.obj)


@cli.command("view")
@click.argument('input_files', required=True, type=click.Path(), nargs=-1)
@click.option("-o", "svs_out", help="Output file, [default: stdout]", required=False, type=click.Path())
@click.option("-f", "--out-format", help="Output format", default="csv", type=click.Choice(["csv", "vcf"]),
              show_default=True)
@click.option("--merge-across", help="Merge records across input samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-within", help="Perform additional merge within input samples, prior to --merge-across",
              default="False", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--separate", help="Keep merged tables separate, adds --post-fix to file names, csv format only",
              default="False", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--post-fix", help="Adds --post-fix to file names, only if --separate is True",
              default="dysgu", type=str, show_default=True)
# @click.option("--filter", help="Pandas compatible filter argument, where d corresponds to output .csv table"
#                                " e.g. 'd['Prob' > 0.5]'",
#               default="", type=str, show_default=False)
# @click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
# @click.option('--exclude', help=".bed file, do not search/call SVs within regions. Overrides include/search",
#               default=None, type=click.Path(exists=True))
# @click.option('--dest', help="Destination folder to use/create for saving results. Defaults to current directory",
#               default=None, type=click.Path())
@click.pass_context
def view_data(ctx, **kwargs):
    """Convert .csv table(s) to .vcf. Merges multiple .csv files into wide .vcf format."""
    # Add arguments to context insert_median, insert_stdev, read_length, out_name
    ctx = apply_ctx(ctx, kwargs)
    return view.view_file(ctx.obj)


@cli.command("test", context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.pass_context
def test_command(ctx, **kwargs):
    """Run dysgu tests"""

    tests_path = os.path.dirname(__file__) + "/tests"

    runner = CliRunner()
    t = [tests_path + '/small.dysgu.sam']
    click.echo(t)
    result = runner.invoke(dysgu_aligner, t)

    t = [tests_path + '/small.bam', '--post-fix', 'dysgu_test']
    click.echo(t)
    result = runner.invoke(get_reads, t)

    t = [tests_path + '/small.bam', '--post-fix', 'dysgu_test']
    click.echo(t)
    result = runner.invoke(find_reads, t)

    t = [tests_path + '/small.bam', '--template-size', '350,100', '-o', './small.dysgu_test.csv']
    click.echo(t)
    result = runner.invoke(call_events, t)

    click.echo("Done", err=True)

