# cython: language_level=3

from __future__ import absolute_import
import click
import os
from sys import argv
import shutil
import time
from multiprocessing import cpu_count
from subprocess import Popen, PIPE
from importlib.metadata import version
import warnings
from dysgu import cluster, view, sv2bam, filter_normals
import datetime
import logging
import shlex
from dysgu.map_set_utils import echo


warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

cpu_range = click.IntRange(min=1, max=cpu_count())


defaults = {
            "clip_length": 15,
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
            "trust_ins_len": "True",
            "symbolic_sv_size": 50_000,
            "spd": 0.3,
            "sd": 0.8,
            "thresholds": "0.45,0.45,0.45,0.45,0.45",
            "divergence": "0.02"
            }


presets = {"nanopore-r9": {"mq": 1,
                        "min_support": "auto",
                        "dist_norm": 900,
                        "max_cov": 150,
                        "pl": "nanopore",
                        "remap": "False",
                        "clip_length": -1,
                        "trust_ins_len": "False",
                        "sd": 0.6,
                        "divergence": "auto",
                        "compression": "wb3",
                        },
            "nanopore-r10": {"mq": 1,
                        "min_support": "auto",
                        "dist_norm": 600,
                        "max_cov": 150,
                        "pl": "nanopore",
                        "remap": "False",
                        "clip_length": -1,
                        "trust_ins_len": "False",
                        "sd": 0.35,
                        "thresholds": "0.35,0.35,0.35,0.35,0.35",
                        "compression": "wb3",
                        },
           "pacbio-sequel2": {"mq": 1,
                      "min_support": "auto",
                      "dist_norm": 600,
                      "max_cov": 150,
                      "pl": "pacbio",
                      "remap": "False",
                      "clip_length": -1,
                      "trust_ins_len": "True",
                      "sd": 0.45,
                      "compression": "wb3",
                      },
           "pacbio-revio": {"mq": 1,
                      "min_support": "auto",
                      "dist_norm": 600,
                      "max_cov": 150,
                      "pl": "pacbio",
                      "remap": "False",
                      "clip_length": -1,
                      "trust_ins_len": "True",
                      "sd": 0.35,
                      "thresholds": "0.25,0.25,0.25,0.25,0.25",
                      "compression": "wb3",
                      },
           "pe": {"mq": defaults["mq"],
                  "min_support": defaults["min_support"],
                  "dist_norm": defaults["dist_norm"],
                  "max_cov": defaults["max_cov"],
                  "pl": defaults["pl"],
                  "remap": defaults["remap"],
                  "trust_ins_len": defaults["trust_ins_len"],
                  },
           }

new_options_set = {}


def add_option_set(ctx, param, value):
    new_options_set[param.name] = value


def show_params():
    logging.info(" ".join(argv[1:]))


def apply_preset(kwargs):
    if kwargs["mode"] == "pacbio": 
        logging.warning("Using --mode pacbio is deprecated. Use 'pacbio-sequel2' or 'pacbio-revio' instead. Mode will be set as 'pacbio-revio'")
        kwargs["mode"] = "pacbio-revio"
    elif kwargs == ["nanopore"]:
        logging.warning("Using --mode nanopore is deprecated. Use 'nanopore-r9' or 'nanopore-r10' instead. Mode will be set as 'nanopore-r10'")
        kwargs["mode"] = "nanopore-r10"
    if kwargs["mode"] != "pe":
        kwargs["paired"] = "False"
    options = new_options_set
    options = {k: v for k, v in options.items() if v is not None}
    for k, v in defaults.items():
        if kwargs.get(k) is None:
            kwargs[k] = v
    kwargs.update(presets[kwargs["mode"]].items())
    kwargs.update(options)


dysgu_version = version("dysgu")

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
        if not temp_dir:
            raise ValueError("Working directory is empty")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
        elif "overwrite" in args and not args["overwrite"]:
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
@click.option("--transcripts", help="A gff3 file of known transcripts. Required for RNAseq data.",
              required=False, type=click.Path(), hidden=True)
@click.option("--sites", help="A vcf file of known variant sites. All sites will be genotyped in the output vcf",
              required=False, type=click.Path())
@click.option("--sites-prob", help="Prior probability that a matching variant in --sites is true",
              required=False, type=click.FloatRange(0, 1), default=0.6, show_default=True)
@click.option('--sites-pass-only', help="Only add variants from sites that have PASS",
              default="True", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option('--ignore-sample-sites', help="If --sites is multi-sample, ignore variants from the input file SV-ALIGNS",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
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
@click.option('--mode', help=f"Type of input reads. Multiple options are set, overrides other options. "
                             f"| pacbio-sequel2: --mq {presets['pacbio-sequel2']['mq']} --paired False --min-support '{presets['pacbio-sequel2']['min_support']}' --max-cov {presets['pacbio-sequel2']['max_cov']} --dist-norm {presets['pacbio-sequel2']['dist_norm']} --trust-ins-len True --sd {presets['pacbio-sequel2']['sd']} --compression wb3. "
                             f"| pacbio-revio: --mq {presets['pacbio-revio']['mq']} --paired False --min-support '{presets['pacbio-revio']['min_support']}' --max-cov {presets['pacbio-revio']['max_cov']} --dist-norm {presets['pacbio-revio']['dist_norm']} --trust-ins-len True --thresholds {presets['pacbio-revio']['thresholds']} --sd {presets['pacbio-revio']['sd']} --compression wb3. "
                             f"| nanopore-r9: --mq {presets['nanopore-r9']['mq']} --paired False --min-support '{presets['nanopore-r9']['min_support']}' --max-cov {presets['nanopore-r9']['max_cov']} --dist-norm {presets['nanopore-r9']['dist_norm']} --trust-ins-len False --sd {presets['nanopore-r9']['sd']} --divergence {presets['nanopore-r9']['divergence']} --compression wb3. "
                             f"| nanopore-r10: --mq {presets['nanopore-r10']['mq']} --paired False --min-support '{presets['nanopore-r10']['min_support']}' --max-cov {presets['nanopore-r10']['max_cov']} --dist-norm {presets['nanopore-r10']['dist_norm']} --trust-ins-len False --thresholds {presets['nanopore-r10']['thresholds']} --sd {presets['nanopore-r10']['sd']} --compression wb3",
              default="pe", type=click.Choice(["pe", "pacbio-sequel2", "pacbio-revio", "nanopore-r9", "nanopore-r10", "pacbio", "nanopore"]), show_default=True)
@click.option('--pl', help=f"Type of input reads  [default: {defaults['pl']}]",
              type=click.Choice(["pe", "pacbio", "nanopore"]), callback=add_option_set)
@click.option('--clip-length', help="Minimum soft-clip length, >= threshold are kept. Set to -1 to ignore [default: {deafults['clip_length']}]", type=int, callback=add_option_set)
@click.option('--max-cov', help=f"Genomic regions with coverage > max-cov discarded."
                                f" Use 'auto' to estimate a value from the alignment index file [default: {defaults['max_cov']}]. Set to -1 to ignore",
              type=str, callback=add_option_set)
@click.option('--max-tlen', help="Maximum template length to consider when calculating paired-end template size",
              default=defaults["max_tlen"], type=int, show_default=True)
@click.option('--min-support', help=f"Minimum number of reads per SV  [default: {defaults['min_support']}]", type=str, callback=add_option_set)
@click.option('--min-size', help="Minimum size of SV to report",
              default=defaults["min_size"], type=int, show_default=True)
@click.option('--mq', help=f"Minimum map quality < threshold are discarded  [default: {defaults['mq']}]",
              type=int, callback=add_option_set)
@click.option('--dist-norm', help=f"Distance normalizer  [default: {defaults['dist_norm']}]", type=float, callback=add_option_set)
@click.option('--spd', help=f"Span position distance [defaults: {defaults['spd']}]", type=float, callback=add_option_set)
@click.option('--sd', help=f"Span distance, only SV span is considered, lower values separate multi-allelic sites [default={defaults['sd']}", type=float, callback=add_option_set)
@click.option('--search-depth', help="Search through this many local reads for matching SVs. Increase this to identify low frequency events", default=20, type=float, show_default=True)
@click.option('--trust-ins-len', help=f"Trust insertion length from cigar, for high error rate reads use False  [default: {defaults['trust_ins_len']}]",
              type=str, callback=add_option_set)
@click.option('--length-extend', help=f"Extend SV length if any nearby gaps found with length >= length-extend. Ignored for paired-end reads", type=int, default=15, show_default=True)
@click.option('--divergence', help=f"Threshold used for ignoring divergent ends of alignments. Ignored for paired-end reads. Use 'auto' to try to infer for noisy reads [default: {defaults['divergence']}]", type=str, callback=add_option_set)
@click.option("-I", "--template-size", help="Manually set insert size, insert stdev, read_length as 'INT,INT,INT'",
              default="", type=str, show_default=False)
@click.option('--search', help=".bed file, limit search to regions", default=None, type=click.Path(exists=True))
@click.option('--exclude', help=".bed file, do not search/call SVs within regions. Takes precedence over --search",
              default=None, type=click.Path(exists=True))
@click.option('--regions', help="bed file of target regions, used for labelling events", default=None, type=click.Path(exists=True))
@click.option('--regions-mm-only', help="If --regions is provided, only use minimizer clustering within --regions. Useful for high coverage targeted sequencing",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("--buffer-size", help="Number of alignments to buffer", default=defaults["buffer_size"],
              type=int, show_default=True)
@click.option("--merge-within", help="Try and merge similar events, recommended for most situations",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--drop-gaps", help="Drop SVs near gaps +/- 250 bp of Ns in reference",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-dist", help="Attempt merging of SVs below this distance threshold. Default for paired-end data is (insert-median + 5*insert_std) for paired reads, or 2000 bp for single-end reads",
              default=None, type=int, show_default=False)
@click.option("--paired", help="Paired-end reads or single", default="True", show_default=True,
              type=click.Choice(["True", "False"]))
@click.option("--contigs", help="Generate consensus contigs for each side of break and use sequence-based metrics in model scoring", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-v', '--verbosity', help="0 = no contigs in output, 1 = output contigs for variants without ALT sequence called, 2 = output all contigs",
              default='1', type=click.Choice(['0', '1', '2']), show_default=True)
@click.option("--diploid", help="Use diploid model for scoring variants. Use 'False' for non-diploid or poly clonal samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--remap", help=f"Try and remap anomalous contigs to find additional small SVs  [default: {defaults['remap']}]", type=str, callback=add_option_set)
@click.option("--no-phase", help="Do not use HP haplotagged reads to phase variants", default=False, is_flag=True, flag_value=True, show_default=True)
@click.option("--metrics", help="Output additional metrics for each SV", default=False, is_flag=True, flag_value=True, show_default=True)
@click.option("--keep-small", help="Keep SVs < min-size found during re-mapping", default=False, is_flag=True, flag_value=True, show_default=False)
@click.option("--symbolic-sv-size", help=f"Use symbolic representation if SV >= this size. Set to -1 to use symbolic-only representation [default={defaults['symbolic_sv_size']}]", type=int, callback=add_option_set)
@click.option("--low-mem", help="Use less memory but more temp disk space", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-c", "--clean", help="Remove temp files and working directory when finished", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--thresholds", help=f"Probability threshold to label as PASS for 'DEL,INS,INV,DUP,TRA' [default: {defaults['thresholds']}", type=str, callback=add_option_set)
@click.pass_context
def run_pipeline(ctx, **kwargs):
    """Run the dysgu pipeline. Important parameters are --mode, --diploid, --min-support, --min-size, --max-cov"""
    # Add arguments to context
    t0 = time.time()
    logging.info("[dysgu-run] Version: {}".format(dysgu_version))

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
    tmp_file_name = f"{dest}/{bname if bname != '-' else os.path.basename(kwargs['working_directory'])}.{pfix}.bam"
    ctx.obj["output"] = tmp_file_name
    ctx.obj["reads"] = "None"

    max_cov_value = sv2bam.process(ctx.obj)
    ctx.obj["max_cov"] = max_cov_value
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
@click.option("--reference", help="Reference file for opening cram files",
              show_default=False, default="", required=False, type=click.Path())
@click.option("-t", "--transcripts", help="A gff3 file of known transcripts. Required for RNAseq data.",
              required=False, type=click.Path(), hidden=True)
@click.option('--pfix', help="Post-fix to add to temp alignment files",
              default="dysgu_reads", type=str)
@click.option("-o", "--output", help="Output reads, discordant, supplementary and soft-clipped reads to file. ",
              type=click.Path(), required=False)
@click.option("--compression", help="Set output bam compression level. Default is uncompressed",
              show_default=True, default="wb0", type=str)
@click.option("-a", "--write_all", help="Write all alignments from SV-read template to temp file", is_flag=True, flag_value=True,
              show_default=True, default=False)
@click.option('--clip-length', help="Minimum soft-clip length, >= threshold are kept. Set to -1 to ignore", type=int, default=15, show_default=True)
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
@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=True, default=False)
@click.option('--pl', help=f"Type of input reads",
              type=click.Choice(["pe", "pacbio", "nanopore"]), default="pe", show_default=True)
@click.pass_context
def get_reads(ctx, **kwargs):
    """Filters input bam/cram for read-pairs that are discordant or have a soft-clip of length > '--clip-length',
    saves bam file in WORKING_DIRECTORY"""
    logging.info("[dysgu-fetch] Version: {}".format(dysgu_version))
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
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option('--ignore-sample-sites', help="If --sites is multi-sample, ignore variants from the input file SV-ALIGNS",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option('--parse-probs', help="Parse INFO:MeanPROB or FORMAT:PROB instead of using --sites-p",
              default="False", type=click.Choice(["True", "False"]),
              show_default=True)
@click.option("--all-sites", help="Output a genotype for all variants in --sites (including homozygous reference 0/0)",
              required=False, default="False", type=click.Choice(["True", "False"]))
@click.option('--pfix', help="Post-fix of temp alignment file (used when a working-directory is provided instead of "
                             "sv-aligns)",
              default="dysgu_reads", type=str, required=False)
@click.option('--mode', help=f"Type of input reads. Multiple options are set, overrides other options. "
                             f"| pacbio-sequel2: --mq {presets['pacbio-sequel2']['mq']} --paired False --min-support '{presets['pacbio-sequel2']['min_support']}' --max-cov {presets['pacbio-sequel2']['max_cov']} --dist-norm {presets['pacbio-sequel2']['dist_norm']} --trust-ins-len True --sd {presets['pacbio-sequel2']['sd']}. "
                             f"| pacbio-revio: --mq {presets['pacbio-revio']['mq']} --paired False --min-support '{presets['pacbio-revio']['min_support']}' --max-cov {presets['pacbio-revio']['max_cov']} --dist-norm {presets['pacbio-revio']['dist_norm']} --trust-ins-len True --thresholds {presets['pacbio-revio']['thresholds']} --sd {presets['pacbio-revio']['sd']}. "
                             f"| nanopore-r9: --mq {presets['nanopore-r9']['mq']} --paired False --min-support '{presets['nanopore-r9']['min_support']}' --max-cov {presets['nanopore-r9']['max_cov']} --dist-norm {presets['nanopore-r9']['dist_norm']} --trust-ins-len False --sd {presets['nanopore-r9']['sd']} --divergence {presets['nanopore-r9']['divergence']}. "
                             f"| nanopore-r10: --mq {presets['nanopore-r10']['mq']} --paired False --min-support '{presets['nanopore-r10']['min_support']}' --max-cov {presets['nanopore-r10']['max_cov']} --dist-norm {presets['nanopore-r10']['dist_norm']} --trust-ins-len False --thresholds {presets['nanopore-r10']['thresholds']} --sd {presets['nanopore-r10']['sd']}",
              default="pe", type=click.Choice(["pe", "pacbio-sequel2", "pacbio-revio", "nanopore-r9", "nanopore-r10", "pacbio", "nanopore"]), show_default=True)
@click.option('--pl', help=f"Type of input reads  [default: {defaults['pl']}]",
              type=click.Choice(["pe", "pacbio", "nanopore"]), callback=add_option_set)
@click.option('--clip-length', help="Minimum soft-clip length, >= threshold are kept. Set to -1 to ignore [default: {deafults['clip_length']}]", type=int, callback=add_option_set)
@click.option('--max-cov', help=f"Regions with > max-cov that do no overlap 'include' are discarded."
                                f" Use 'auto' to estimate a value from the alignment index file [default: {defaults['max_cov']}]. Regions with > max-cov that do no overlap 'include' are discarded. Set to -1 to ignore.",
              type=str, callback=add_option_set)
@click.option('--max-tlen', help="Maximum template length to consider when calculating paired-end template size",
              default=defaults["max_tlen"], type=int, show_default=True)
@click.option('--min-support', help=f"Minimum number of reads per SV  [default: {defaults['min_support']}]", type=str, callback=add_option_set)
@click.option('--min-size', help="Minimum size of SV to report",
              default=defaults["min_size"], type=int, show_default=True)
@click.option('--mq', help=f"Minimum map quality < threshold are discarded  [default: {defaults['mq']}]",
              type=int, callback=add_option_set)
@click.option('--dist-norm', help=f"Distance normalizer  [default: {defaults['dist_norm']}]", type=float, callback=add_option_set)
@click.option('--spd', help=f"Span position distance [defaults: {defaults['spd']}]", type=float, callback=add_option_set)

@click.option('--sd', help=f"Span distance, only SV span is considered, lower values separate multi-allelic sites [default={defaults['sd']}", type=float, callback=add_option_set)
@click.option('--search-depth', help="Search through this many local reads for matching SVs. Increase this to identify low frequency events", default=20, type=float, show_default=True)
@click.option('--trust-ins-len', help=f"Trust insertion length from cigar, for high error rate reads use False  [default: {defaults['trust_ins_len']}]", type=str, callback=add_option_set)
@click.option('--length-extend', help=f"Extend SV length if any nearby gaps found with length >= length-extend. Ignored for paired-end reads", type=int, default=15, show_default=True)
@click.option('--divergence', help=f"Threshold used for ignoring divergent ends of alignments. Ignored for paired-end reads. Use 'auto' to try to infer for noisy reads [default: {defaults['divergence']}]", type=str, callback=add_option_set)
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
                                   "reads, or 2000 bp for single-end reads",
              default=None, type=int, show_default=False)
@click.option("--paired", help="Paired-end reads or single", default="True", show_default=True,
              type=click.Choice(["True", "False"]))
@click.option("--contigs", help="Generate consensus contigs for each side of break and use sequence-based metrics in model scoring", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-v', '--verbosity', help="0 = no contigs in output, 1 = output contigs for variants without ALT sequence called, 2 = output all contigs",
              default='1', type=click.Choice(['0', '1', '2']), show_default=True)
@click.option("--diploid", help="Use diploid model for scoring variants. Use 'False' for non-diploid or poly clonal samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--remap", help=f"Try and remap anomalous contigs to find additional small SVs  [default: {defaults['remap']}]", type=str, callback=add_option_set)
@click.option("--no-phase", help="Do not use HP haplotagged reads to phase variants", default=False, is_flag=True, flag_value=True, show_default=True)
@click.option("--metrics", help="Output additional metrics for each SV", default=False, is_flag=True, flag_value=True, show_default=True)
@click.option("--keep-small", help="Keep SVs < min-size found during re-mapping", default=False, is_flag=True, flag_value=True, show_default=False)
@click.option("--symbolic-sv-size", help=f"Use symbolic representation if SV >= this size. Set to -1 to use symbolic-only representation [default={defaults['symbolic_sv_size']}]", type=int, callback=add_option_set)
@click.option("--low-mem", help="Use less memory but more temp disk space", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-x", "--overwrite", help="Overwrite temp files", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("-c", "--clean", help="Remove temp files and working directory when finished", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--thresholds", help=f"Probability threshold to label as PASS for 'DEL,INS,INV,DUP,TRA' [default: {defaults['thresholds']}", type=str, callback=add_option_set)
@click.pass_context
def call_events(ctx, **kwargs):
    """Call structural variants from bam alignment file/stdin"""
    logging.info("[dysgu-call] Version: {}".format(dysgu_version))
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
    if kwargs["sv_aligns"].endswith(".cram"):
        raise ValueError("Cram files are not supported for 'call' command, please use 'run' pipeline instead")
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
@click.argument('input_files', required=False, type=click.Path(), nargs=-1)
@click.option("-i", "--input-list", help="Input list of file paths, one line per file ", required=False, type=click.Path(exists=True))
@click.option("-o", "--svs-out", help="Output file, [default: stdout]", required=False, type=click.Path())
@click.option("-f", "--out-format", help="Output format", default="vcf", type=click.Choice(["csv", "vcf"]),
              show_default=True)
@click.option("-p", "--procs", help="Number of processors to use when merging, requires --wd option to be supplied", type=cpu_range, default=1, show_default=True)
@click.option("-d", "--wd", help="Working directory to use/create when merging", type=click.Path(exists=False), required=False)
@click.option("-c", "--clean", help="Remove working directory when finished", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--progress", help="Prints detailed progress information",  is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--cohort-update", help="Updated this cohort file with new calls from input_files", required=False, type=click.Path(exists=True))
@click.option("--collapse-nearby", help="Merges more aggressively by collapsing nearby SV",
              default="True", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-across", help="Merge records across input samples", default="True",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-method", help="Method of merging using --merge-across. Progressive is suitable for large cohorts. Auto will switch to progressive for >4 samples",
              default="auto", type=click.Choice(["auto", "all-vs-all", "progressive"]), show_default=True)
@click.option("--merge-within", help="Perform additional merge within input samples, prior to --merge-across",
              default="False", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--merge-dist", help="Distance threshold for merging",
              default=500, type=int, show_default=True)
@click.option("--max-comparisons", help="Compare each event with up to --max-comparisons local SVs",
              default=20, type=int, show_default=True)
@click.option("--separate", help="Keep merged tables separate, adds --post-fix to file names, csv format only",
              default="False", type=click.Choice(["True", "False"]), show_default=True)
@click.option("--post-fix", help="Adds --post-fix to file names, only if --separate is True",
              default="dysgu", type=str, show_default=True)
@click.option("--add-kind", help="Add region-overlap 'kind' to vcf output", default="False",
              type=click.Choice(["True", "False"]), show_default=True)
@click.option('-v', '--verbosity', help="0 = no contigs in output, 1 = output contigs for variants without ALT sequence called, 2 = output all contigs",
              default='1', type=click.Choice(['0', '1', '2']), show_default=True)
@click.pass_context
def view_data(ctx, **kwargs):
    """Merge vcf/csv variant files"""
    logging.info("[dysgu-merge] Version: {}".format(dysgu_version))
    ctx = apply_ctx(ctx, kwargs)
    if ctx.obj["wd"]:
        make_wd(ctx.obj)
    if not ctx.obj["input_files"] and not ctx.obj["input_list"]:
        logging.error("Please supply either INPUT_FILES as an argument or use the --input-list option")
        quit()
    view.view_file(ctx.obj)


@cli.command("filter")
@click.argument('input_vcf', required=True, type=click.Path())
@click.argument('normal_bams', required=False, type=click.Path(), nargs=-1)
@click.option('--reference', help="Reference for cram input files", required=False, type=click.Path(exists=True))
@click.option("-o", "--svs-out", help="Output file, [default: stdout]", required=False, type=click.Path())
@click.option("-n", "--normal-vcf", help="Vcf file for normal sample, or panel of normals. The SM tag of input bams is used to ignore the input_vcf for multi-sample vcfs", required=False, type=click.Path(), multiple=True)
@click.option("-p", "--procs", help="Reading threads for normal_bams", type=cpu_range, default=1, show_default=True)
@click.option("-f", "--support-fraction", help="Minimum threshold support fraction / coverage (SU/COV)", type=float, default=0.1, show_default=True)
@click.option("--target-sample", help="If input_vcf if multi-sample, use target-sample as input", required=False, type=str, default="")
@click.option("--keep-all", help="All SVs classified as normal will be kept in the output, labelled as filter=normal", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--ignore-read-groups", help="Ignore ReadGroup RG tags when parsing sample names. Filenames will be used instead", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--max-divergence", help="Remove SV if normal_bam displays divergence > max-divergence at same location", default=0.1, type=float, show_default=True)
@click.option("--min-prob", help="Remove SV with PROB value < min-prob", default=0.1, type=float, show_default=True)
@click.option("--min-mapq", help="Remove SV with mean mapqq < min-mapq", default=10, type=float, show_default=True)
@click.option("--pass-prob", help="Re-label SV as PASS if PROB value >= pass-prob", default=1.0, type=float, show_default=True)
@click.option("--interval-size", help="Interval size for searching normal-vcf/normal-bams", default=1000, type=int, show_default=True)
@click.option("--random-bam-sample", help="Choose N random normal-bams to search. Use -1 to ignore", default=-1, type=int, show_default=True)
@click.pass_context
def filter_normal(ctx, **kwargs):
    """Filter a vcf generated by dysgu.
    Unique SVs can be found in the input_vcf by supplying a --normal-vcf (single or multi-sample),
    and normal bam files. Bam/vcf samples with the same name as the input_vcf will be ignored"""
    logging.info("[dysgu-filter] Version: {}".format(dysgu_version))
    ctx = apply_ctx(ctx, kwargs)
    filter_normals.run_filtering(ctx.obj)


@cli.command("test")
@click.option("--verbose", help="Show output of test commands",
              is_flag=True, flag_value=True, show_default=False, default=False)
@click.pass_context
def test_command(ctx, **kwargs):
    """Run dysgu tests"""
    pwd = os.getcwd()
    logging.info("[dysgu-test] Version: {}".format(dysgu_version))
    if kwargs["verbose"]:
        logging.info(f"Current directory: {pwd}") 
    tests_path = os.path.dirname(__file__) + "/tests"
    tests = list()

    tests.append(["dysgu fetch",
                  "-x ",
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False",
                  "-o " + pwd + '/test.dysgu{}.vcf'.format(dysgu_version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False",
                  "--regions " + tests_path + '/targets.bed',
                  "-o " + pwd + '/test_regions.dysgu{}.vcf'.format(dysgu_version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test2',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False --procs 2",
                  "--regions " + tests_path + '/targets.bed',
                  "-o " + pwd + '/test_regions.dysgu{}.vcf'.format(dysgu_version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test2',
                  tests_path + '/small.bam'])
    tests.append(["dysgu run",
                  "-x --drop-gaps False --mode pacbio-sequel2",
                  "-o " + pwd + '/test2.dysgu{}.vcf'.format(dysgu_version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu call",
                  "-x --drop-gaps False",
                  "-o " + pwd + '/test2.dysgu{}.vcf'.format(dysgu_version),
                  tests_path + '/ref.fa',
                  pwd + '/wd_test',
                  tests_path + '/small.bam'])
    tests.append(["dysgu merge",
                  pwd + '/test.dysgu{}.vcf'.format(dysgu_version),
                  pwd + '/test2.dysgu{}.vcf'.format(dysgu_version),
                  "-o " + pwd + '/test.merge.dysgu{}.vcf'.format(dysgu_version)])
    for t in tests:
        c = " ".join(t)
        process = Popen(shlex.split(c), stdout=PIPE, stderr=PIPE, text=True)
        if kwargs["verbose"]:
            for line in process.stderr:
                click.echo(line.strip(), err=True)
        process.wait()
        if process.returncode != 0:
            logging.warning(f"WARNING: Command failed with return code {process.returncode}\n{c}")
        else:
            click.echo("PASS: " + c + "\n", err=True)

    logging.info("Run test complete")


if __name__ == "__main__":
    cli()
