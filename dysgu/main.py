from __future__ import absolute_import
import click
import os
from sys import argv
import shutil
import time
from multiprocessing import cpu_count
from subprocess import run
import pkg_resources
import warnings
from dysgu import cluster, view, sv2bam
import datetime
import logging
from dysgu.map_set_utils import echo


warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

cpu_range = click.IntRange(min=1, max=cpu_count())

defaults = {
            "clip_length": 15,
            "output": "-",
            "svs_out": "-",
            "max_cov": 200,
            "buffer_size": 0,
            "min_support": 3,
            "min_size": 30,
            "model": None,
            "max_tlen": 1000,
            "z_depth": 2,
            "z_breadth": 2,
            "dist_norm": 100,
            "mq": 1,
            "regions_only": "False",
            "soft_search": "True",
            "default_long_support": 2,
            "pl": "pe",
            "remap": "True",
            "drop_gaps": "True",
            "trust_ins_len": "True"
            }


presets = {"nanopore": {"mq": 20,
                        "min_support": 2,
                        "dist_norm": 900,
                        "max_cov": 150,
                        "pl": "nanopore",
                        "remap": "False",
                        "clip_length": -1,
                        "trust_ins_len": "False"
                        },
           "pacbio": {"mq": 20,
                      "min_support": 2,
                      "dist_norm": 600,
                      "max_cov": 150,
                      "pl": "pacbio",
                      "remap": "False",
                      "clip_length": -1,
                      "trust_ins_len": "True"
                      },
           "pe": {"mq": defaults["mq"],
                  "min_support": defaults["min_support"],
                  "dist_norm": defaults["dist_norm"],
                  "max_cov": defaults["max_cov"],
                  "pl": defaults["pl"],
                  "remap": defaults["remap"],
                  "trust_ins_len": defaults["trust_ins_len"]
                  },
           }

new_options_set = {}


def add_option_set(ctx, param, value):
    new_options_set[param.name] = value


def show_params():
    logging.info(" ".join(argv[1:]))


def apply_preset(kwargs):

    if kwargs["mode"] == "nanopore" or kwargs["mode"] == "pacbio":
        kwargs["paired"] = "False"

    for k, v in presets[kwargs["mode"]].items():
        if k in new_options_set and new_options_set[k] is not None:

            # use the presents option if new_options_set is default
            # if new_options_set[k] == defaults[k]:
            #     kwargs[k] = v

            kwargs[k] = new_options_set[k]
        else:
            kwargs[k] = v

    for k, v in defaults.items():
        if k in kwargs and kwargs[k] is None:
            kwargs[k] = v


version = pkg_resources.require("dysgu")[0].version


logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)


def apply_ctx(ctx, kwargs):
    ctx.ensure_object(dict)
    if len(ctx.obj) == 0:  # When invoked from cmd line, else run was invoked from test function
        for k, v in list(defaults.items()) + list(kwargs.items()):
            ctx.obj[k] = v
    return ctx


def make_wd(args, call_func=False):
    if "working_directory" in args:
        temp_dir = args["working_directory"]

        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
        elif not args["overwrite"]:
            if (call_func and args["ibam"] is None) or not call_func:
                raise ValueError("Working directory already exists. Add -x / --overwrite=True to proceed, "
                                 "or supply --ibam to re-use temp files in working directory")

    else:
        temp_dir = args["wd"]
        if temp_dir is None:
            return -1
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
        elif not args["overwrite"]:
            raise ValueError("Working directory already exists. Add -x / --overwrite=True to proceed")


@click.group(chain=False, invoke_without_command=False)
@click.version_option()
def cli():
    """Dysgu-SV is a set of tools calling structural variants from bam/cram files"""
    pass


@cli.command("run")
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('working_directory', required=True, type=click.Path())
@click.argument('bam', required=True, type=click.Path(exists=False))
@click.option("--sites", help="A vcf file of known variant sites. All sites will be genotyped in the output vcf",
              required=False, type=click.Path())
@click.option("--sites-prob", help="Prior probability that a matching variant in --sites is true",
              required=False, type=click.FloatRange(0, 1), default=0.6, show_default=True)
@click.option('--sites-pass-only', help="Only add variants from sites that have PASS",
              default="True", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option('--parse-probs', help="Parse INFO:MeanPROB or FORMAT:PROB instead of using --sites-p",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("--all-sites", help="Output a genotype for all variants in --sites (including homozygous reference 0/0)",
              required=False, default="False", type=click.Choice(["True", "False"]))
@click.option('--pfix', help="Post-fix to add to temp alignment files", default="dysgu_reads", type=str)
@click.option("-o", "--svs-out", help="Output file, [default: stdout]", required=False, type=click.Path())
@click.option("-f", "--out-format", help="Output format", default="vcf", type=click.Choice(["csv", "vcf"]),
              show_default=True)
@click.option("-a", "--write_all", help="Write all alignments from SV-read template to temp file", is_flag=True, flag_value=True,
              show_default=True, default=False)
@click.option("--compression", help="Set temp file bam compression level. Default is uncompressed",
              show_default=True, default="wb0", type=str)
@click.option("-p", "--procs", help="Number of cpu cores to use", type=cpu_range, default=1,
              show_default=True)
@click.option('--mode', help="Type of input reads. Multiple options are set, overrides other options. "
                             "pacbio: --mq 20 --paired False --min-support 2 --max-cov 150 --dist-norm 200 --trust-ins-len True. "
                             "nanopore: --mq 20 --paired False --min-support 2 --max-cov 150 --dist-norm 900 --trust-ins-len False",
              default="pe", type=click.Choice(["pe", "pacbio", "nanopore"]), show_default=True)
@click.option('--pl', help=f"Type of input reads  [default: {defaults['pl']}]",
              type=click.Choice(["pe", "pacbio", "nanopore"]), callback=add_option_set)
@click.option('--clip-length', help="Minimum soft-clip length, >= threshold are kept. Set to -1 to ignore", default=defaults["clip_length"],
              type=int, show_default=True)
@click.option('--max-cov', help=f"Genomic regions with coverage > max-cov discarded."
                                f" Use 'auto' to estimate a value from the alignment index file [default: {defaults['max_cov']}]. Set to -1 to ignore",
              type=str, callback=add_option_set)
@click.option('--max-tlen', help="Maximum template length to consider when calculating paired-end template size",
              default=defaults["max_tlen"], type=int, show_default=True)
@click.option('--min-support', help=f"Minimum number of reads per SV  [default: {defaults['min_support']}]", type=int, callback=add_option_set)
@click.option('--min-size', help="Minimum size of SV to report",
              default=defaults["min_size"], type=int, show_default=True)
@click.option('--mq', help=f"Minimum map quality < threshold are discarded  [default: {defaults['mq']}]",
              type=int, callback=add_option_set)
@click.option('--dist-norm', help=f"Distance normalizer  [default: {defaults['dist_norm']}]", type=float, callback=add_option_set)
@click.option('--spd', help="Span position distance", default=0.3, type=float, show_default=True)
@click.option('--trust-ins-len', help=f"Trust insertion length from cigar, for high error rate reads use False  [default: {defaults['trust_ins_len']}]",
              type=str, callback=add_option_set)
@click.option("-I", "--template-size", help="Manually set insert size, insert stdev, read_length as 'INT,INT,INT'",
              default="", type=str, show_default=False)
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Takes precedence over --search",
              default=None, type=click.Path(exists=True))
@click.option('--regions', help="bed file of target regions, used for labelling events", default=None, type=click.Path(exists=True))
@click.option('--regions-only', help="If --regions is provided, call only events within target regions",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option('--regions-mm-only', help="If --regions is provided, only use minimizer clustering within --regions. Useful for high coverage targeted sequencing",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("--buffer-size", help="Number of alignments to buffer", default=defaults["buffer_size"],
              type=int, show_default=True)
@click.option("--merge-within", help="Try and merge similar events, recommended for most situations",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--drop-gaps", help="Drop SVs near gaps +/- 250 bp of Ns in reference",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-dist", help="Attempt merging of SVs below this distance threshold. Default for paired-end data is (insert-median + 5*insert_std) for paired"
                                   "reads, or 500 bp for single-end reads",
              default=None, type=int, show_default=False)
@click.option("--paired", help="Paired-end reads or single", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--contigs", help="Generate consensus contigs for each side of break and use sequence-based metrics in model scoring", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-v', '--verbosity', help="0 = no contigs in output, 1 = output contigs for variants without ALT sequence called, 2 = output all contigs",
              default='1', type=click.Choice(['0', '1', '2']), show_default=True)
@click.option("--diploid", help="Use diploid model for scoring variants. Use 'False' for non-diploid or poly clonal samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--remap", help=f"Try and remap anomalous contigs to find additional small SVs  [default: {defaults['remap']}]", type=str, callback=add_option_set)
@click.option("--metrics", help="Output additional metrics for each SV", default=False, is_flag=True, flag_value=True, show_default=True)
@click.option("--no-gt", help="Skip adding genotype to SVs", is_flag=True, flag_value=False, show_default=False, default=True)
@click.option("--keep-small", help="Keep SVs < min-size found during re-mapping", default=False, is_flag=True, flag_value=True, show_default=False)
@click.option("--low-mem", help="Use less memory but more temp disk space", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-c", "--clean", help="Remove temp files when finished", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--thresholds", help="Probability threshold to label as PASS for 'DEL,INS,INV,DUP,TRA'", default="0.45,0.45,0.45,0.45,0.45",
              type=str, show_default=True)
@click.pass_context
def run_pipeline(ctx, **kwargs):
    """Run the dysgu pipeline. Important parameters are --mode, --diploid, --min-support, --min-size, --max-cov"""
    # Add arguments to context
    t0 = time.time()
    logging.info("[dysgu-run] Version: {}".format(version))
    make_wd(kwargs)

    apply_preset(kwargs)
    show_params()

    ctx = apply_ctx(ctx, kwargs)

    if kwargs["diploid"] == "False" and kwargs["contigs"] == "False":
        raise ValueError("Only dip=False or contigs=False are supported, not both")

    pfix = kwargs["pfix"]

    dest = os.path.expanduser(kwargs["working_directory"])
    logging.info(f"Destination: {dest}")
    bname = os.path.splitext(os.path.basename(kwargs["bam"]))[0]
    tmp_file_name = f"{dest}/{bname if bname != '-' else dest}.{pfix}.bam"

    # Get SV reads
    if kwargs["regions"] is not None and kwargs["regions_only"] == "True":
        ctx.obj["search"] = kwargs["regions"]

    ctx.obj["output"] = tmp_file_name
    ctx.obj["reads"] = "None"

    max_cov_value = sv2bam.process(ctx.obj)

    ctx.obj["max_cov"] = max_cov_value
    # Call SVs
    if kwargs["bam"] != "-":
        ctx.obj["ibam"] = kwargs["bam"]
    else:
        ctx.obj["ibam"] = None
    ctx.obj["sv_aligns"] = tmp_file_name
    logging.info("Input file is: {}".format(tmp_file_name))
    cluster.cluster_reads(ctx.obj)

    if kwargs["clean"]:
        shutil.rmtree(kwargs["working_directory"])

    logging.info("dysgu run {} complete, time={} h:m:s".format(kwargs["bam"], str(datetime.timedelta(
        seconds=int(time.time() - t0)))))


@cli.command("fetch")
@click.argument('working_directory', required=True, type=click.Path())
@click.argument('bam', required=True, type=click.Path(exists=False))
# @click.option("-f", "--out-format", help="Output format. 'bam' output maintains sort order, "
#                                          "'fq' output is collated by name",
#               default="bam", type=click.Choice(["bam", "fq", "fasta"]),
#               show_default=True)
@click.option("--reference", help="Reference file for opening cram files",
              show_default=False, default="", required=False, type=click.Path())
@click.option('--pfix', help="Post-fix to add to temp alignment files",
              default="dysgu_reads", type=str)
@click.option("-o", "--output", help="Output reads, discordant, supplementary and soft-clipped reads to file. ",
              type=str)
@click.option("--compression", help="Set output bam compression level. Default is uncompressed",
              show_default=True, default="wb0", type=str)
@click.option("-a", "--write_all", help="Write all alignments from SV-read template to temp file", is_flag=True, flag_value=True,
              show_default=True, default=False)
@click.option('--clip-length', help="Minimum soft-clip length, >= threshold are kept. Set to -1 to ignore",
              default=defaults["clip_length"],
              type=int, show_default=True)
@click.option('--mq', help="Minimum map quality < threshold are discarded", default=1,
              type=int, show_default=True)
@click.option('--min-size', help="Minimum size of SV to report",
              default=defaults["min_size"], type=int, show_default=True)
@click.option('--max-cov', help="Genomic regions with coverage > max-cov are discarded. Set to -1 to ignore.",
              default=defaults["max_cov"], type=float, show_default=True)
@click.option("-p", "--procs", help="Compression threads to use for writing bam", type=cpu_range, default=1,
              show_default=True)
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Takes precedence over --search",
              default=None, type=click.Path(exists=True))
# @click.option('--regions', help="bed file of target regions. --max-cov is ignored within these regions", default=None,
#               type=click.Path(exists=True))
@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=True, default=False)
@click.option('--pl', help=f"Type of input reads  [default: {defaults['pl']}]",
              type=click.Choice(["pe", "pacbio", "nanopore"]), callback=add_option_set)
@click.pass_context
def get_reads(ctx, **kwargs):
    """Filters input .bam/.cram for read-pairs that are discordant or have a soft-clip of length > '--clip-length',
    saves bam file in WORKING_DIRECTORY"""
    logging.info("[dysgu-fetch] Version: {}".format(version))
    make_wd(kwargs)

    ctx = apply_ctx(ctx, kwargs)
    return sv2bam.process(ctx.obj)


@cli.command("call")
@click.argument('reference', required=True, type=click.Path(exists=True))
@click.argument('working_directory', required=True, type=click.Path())
@click.argument('sv-aligns', required=False, type=click.Path(exists=False))
@click.option("-b", "--ibam", help="Original input file usef with 'fetch' command, used for calculating insert size parameters",
              show_default=True, default=None, required=False, type=click.Path())
@click.option("-o", "--svs-out", help="Output file [default: stdout]", required=False, type=click.Path())
@click.option("-f", "--out-format", help="Output format", default="vcf", type=click.Choice(["csv", "vcf"]),
              show_default=True)
@click.option("--sites", help="A vcf file of known variant sites. Matching output variants are labelled with 'PASS' plus the ID from --sites",
              required=False, type=click.Path())
@click.option("--sites-prob", help="Prior probability that a matching variant in --sites is true",
              required=False, type=click.FloatRange(0, 1), default=0.6, show_default=True)
@click.option('--sites-pass-only', help="Only add variants from sites that have PASS",
              default="True", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option('--parse-probs', help="Parse INFO:MeanPROB or FORMAT:PROB instead of using --sites-p",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("--all-sites", help="Output a genotype for all variants in --sites (including homozygous reference 0/0)",
              required=False, default="False", type=click.Choice(["True", "False"]))
@click.option('--pfix', help="Post-fix of temp alignment file (used when a working-directory is provided instead of "
                             "sv-aligns)",
              default="dysgu_reads", type=str, required=False)
@click.option('--mode', help="Type of input reads. Multiple options are set, overrides other options"
                             "pacbio/nanopore: --mq 20 --paired False --min-support 2 --max-cov 150", default="pe",
              type=click.Choice(["pe", "pacbio", "nanopore"]), show_default=True)
@click.option('--pl', help=f"Type of input reads  [default: {defaults['pl']}]",
              type=click.Choice(["pe", "pacbio", "nanopore"]), callback=add_option_set)
@click.option('--clip-length', help="Minimum soft-clip length, >= threshold are kept. Set to -1 to ignore", default=defaults["clip_length"],
              type=int, show_default=True)
@click.option('--max-cov', help=f"Regions with > max-cov that do no overlap 'include' are discarded."
                                f" Use 'auto' to estimate a value from the alignment index file [default: {defaults['max_cov']}]. Regions with > max-cov that do no overlap 'include' are discarded. Set to -1 to ignore.",
              type=str, callback=add_option_set)
@click.option('--max-tlen', help="Maximum template length to consider when calculating paired-end template size",
              default=defaults["max_tlen"], type=int, show_default=True)
@click.option('--min-support', help=f"Minimum number of reads per SV  [default: {defaults['min_support']}]", type=int, callback=add_option_set)
@click.option('--min-size', help="Minimum size of SV to report",
              default=defaults["min_size"], type=int, show_default=True)
@click.option('--mq', help=f"Minimum map quality < threshold are discarded  [default: {defaults['mq']}]",
              type=int, callback=add_option_set)
@click.option('--dist-norm', help=f"Distance normalizer  [default: {defaults['dist_norm']}]", type=float, callback=add_option_set)
@click.option('--spd', help="Span position distance", default=0.3, type=float, show_default=True)
@click.option('--trust-ins-len', help=f"Trust insertion length from cigar, for high error rate reads use False  [default: {defaults['trust_ins_len']}]",
              type=str, callback=add_option_set)
@click.option("-I", "--template-size", help="Manually set insert size, insert stdev, read_length as 'INT,INT,INT'",
              default="", type=str, show_default=False)
@click.option('--regions', help="bed file of target regions, used for labelling events", default=None, type=click.Path(exists=True))
@click.option('--regions-only', help="If --regions is provided, call only events within target regions",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option('--regions-mm-only', help="If --regions is provided, only use minimizer clustering within --regions. Useful for high coverage targeted sequencing",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("-p", "--procs", help="Processors to use", type=cpu_range, default=1, show_default=True)
@click.option("--buffer-size", help="Number of alignments to buffer", default=defaults["buffer_size"],
              type=int, show_default=True)
@click.option("--merge-within", help="Try and merge similar events, recommended for most situations",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--drop-gaps", help="Drop SVs near gaps +/- 250 bp of Ns in reference",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-dist", help="Attempt merging of SVs below this distance threshold, default is (insert-median + 5*insert_std) for paired"
                                   "reads, or 500 bp for single-end reads",
              default=None, type=int, show_default=False)
@click.option("--paired", help="Paired-end reads or single", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--contigs", help="Generate consensus contigs for each side of break and use sequence-based metrics in model scoring", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-v', '--verbosity', help="0 = no contigs in output, 1 = output contigs for variants without ALT sequence called, 2 = output all contigs",
              default='1', type=click.Choice(['0', '1', '2']), show_default=True)
@click.option("--diploid", help="Use diploid model for scoring variants. Use 'False' for non-diploid or poly clonal samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--remap", help=f"Try and remap anomalous contigs to find additional small SVs  [default: {defaults['remap']}]", type=str, callback=add_option_set)
@click.option("--metrics", help="Output additional metrics for each SV", default=False, is_flag=True, flag_value=True, show_default=True)
@click.option("--no-gt", help="Skip adding genotype to SVs", is_flag=True, flag_value=False, show_default=False, default=True)
@click.option("--keep-small", help="Keep SVs < min-size found during re-mapping", default=False, is_flag=True, flag_value=True, show_default=False)
@click.option("--low-mem", help="Use less memory but more temp disk space", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-c", "--clean", help="Remove temp files when finished", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--thresholds", help="Probability threshold to label as PASS for 'DEL,INS,INV,DUP,TRA'", default="0.45,0.45,0.45,0.45,0.45",
              type=str, show_default=True)
@click.pass_context
def call_events(ctx, **kwargs):
    """Call structural variants from alignment file/stdin"""

    logging.info("[dysgu-call] Version: {}".format(version))
    make_wd(kwargs, call_func=True)
    if kwargs["sv_aligns"] is None:
        # Try and open from working director
        if os.path.exists(kwargs["working_directory"]):
            bname = os.path.basename(kwargs["working_directory"])
            pth = "{}/{}.{}.bam".format(kwargs["working_directory"], bname, kwargs["pfix"])
            if os.path.exists(pth):
                kwargs["sv_aligns"] = pth
            else:
                raise ValueError("Could not find '{}'".format(kwargs["sv_aligns"]))

    if kwargs["diploid"] == "False" and kwargs["contigs"] == "False":
        raise ValueError("Only dip=False or contigs=False are supported, not both")
    logging.info("Input file is: {}".format(kwargs["sv_aligns"]))

    apply_preset(kwargs)

    if kwargs["max_cov"] == "-1":
        kwargs["max_cov"] = 1e6
    show_params()
    ctx = apply_ctx(ctx, kwargs)

    cluster.cluster_reads(ctx.obj)

    if kwargs["clean"]:
        shutil.rmtree(kwargs["working_directory"])


@cli.command("merge")
@click.argument('input_files', required=True, type=click.Path(), nargs=-1)
@click.option("-o", "svs_out", help="Output file, [default: stdout]", required=False, type=click.Path())
#@click.option('--pon', help="'Pool-of-normals' comma-separated list of indexed alignment files. If query SVs found in any of --pon, SV is removed", default=None, type=str)
#@click.option('--pon-rt', help="Read-types of the --pon samples (comma-separated list). Either pe/pacbio/nanopore", default=None, type=str)
#@click.option('--ref', help="Reference file, required if --pon used", required=False, type=click.Path(exists=True))
#@click.option('--wd', help="Working directory, required if --pon used", required=False, type=click.Path())
@click.option("-f", "--out-format", help="Output format", default="vcf", type=click.Choice(["csv", "vcf"]),
              show_default=True)
@click.option("--merge-across", help="Merge records across input samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-within", help="Perform additional merge within input samples, prior to --merge-across",
              default="False", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-dist", help="Distance threshold for merging",
              default=500, type=int, show_default=True)
@click.option("--separate", help="Keep merged tables separate, adds --post-fix to file names, csv format only",
              default="False", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--post-fix", help="Adds --post-fix to file names, only if --separate is True",
              default="dysgu", type=str, show_default=True)
@click.option("--no-chr", help="Remove 'chr' from chromosome names in vcf output", default="False",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--no-contigs", help="Remove contig sequences from vcf output", default="False",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--add-kind", help="Add region-overlap 'kind' to vcf output", default="False",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-v', '--verbosity', help="0 = no contigs in output, 1 = output contigs for variants without ALT sequence called, 2 = output all contigs",
              default='1', type=click.Choice(['0', '1', '2']), show_default=True)
#@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=False, default=False)
#@click.option("-c", "--clean", help="Remove temp files when finished", is_flag=True, flag_value=True, show_default=False, default=False)
@click.pass_context
def view_data(ctx, **kwargs):
    """Merge .vcf/csv variant files"""
    # Add arguments to context insert_median, insert_stdev, read_length, out_name
    logging.info("[dysgu-merge] Version: {}".format(version))
    # make_wd(kwargs, call_func=True)
    ctx = apply_ctx(ctx, kwargs)
    view.view_file(ctx.obj)
    #if kwargs["wd"] is not None and kwargs["clean"]:
    #    shutil.rmtree(kwargs["wd"])


@cli.command("test")
@click.pass_context
def test_command(ctx, **kwargs):
    """Run dysgu tests"""
    pwd = os.getcwd()
    logging.info("[dysgu-test] Version: {}".format(version))
    tests_path = os.path.dirname(__file__) + "/tests"

    tests = list()
    tests.append(["dysgu fetch",
                  "-x ",
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False",
                  "-o " + pwd + '/test.dysgu{}.vcf'.format(version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False",
                  "--regions " + tests_path + '/targets.bed',
                  "-o " + pwd + '/test_regions.dysgu{}.vcf'.format(version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test2',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False --procs 2",
                  "--regions " + tests_path + '/targets.bed',
                  "-o " + pwd + '/test_regions.dysgu{}.vcf'.format(version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test2',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False --mode pacbio",
                  "-o " + pwd + '/test2.dysgu{}.vcf'.format(version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu call",
                  "-x --drop-gaps False",
                  "-o " + pwd + '/test2.dysgu{}.vcf'.format(version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu merge",
                  pwd + '/test.dysgu{}.vcf'.format(version),
                  pwd + '/test2.dysgu{}.vcf'.format(version),
                  "-o " + pwd + '/test.merge.dysgu{}.vcf'.format(version)])

    for t in tests:
        c = " ".join(t)
        v = run(c, shell=True, capture_output=True, check=True)
        if v.returncode != 0:
            raise RuntimeError(t, "finished with non zero: {}".format(v))
        else:
            click.echo("PASS: " + c, err=True)

    logging.info("Run test complete")
