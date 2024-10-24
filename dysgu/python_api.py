import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

from dysgu import post_call as post_call_metrics
from dysgu.cluster import pipe1
from dysgu.io_funcs import to_vcf
from dysgu.map_set_utils import to_dict, merge_intervals, load_bed
from dysgu.merge_svs import merge_events
from dysgu.view import dotdict, set_numeric, vcf_to_df


def dysgu_default_args():
    """
    Returns the default arguments used by dysgu
    :return: A dict of available arguments
    :rtype: dict
    """
    args = {
                'add_gt': True,
                'add_kind': True,
                'all_sites': False,
                'buffer_size': 10000,
                'clean': False,
                'clip_length': 15,
                'contigs': True,
                'diploid': True,
                'dist_norm': 100,
                'divergence': 0.02,
                'drop_gaps': True,
                'ibam': None,
                'ignore_sample_sites': True,
                'keep_small': False,
                'length_extend': 15,
                'low_mem': False,
                'max_cov': 200,
                'max_tlen': 1000,
                'merge_dist': None,
                'merge_within': True,
                'metrics': True,
                'min_size': 30,
                'min_support': 3,
                'mode': 'pe',
                'model': None,
                'mq': 1,
                'overwrite': True,
                'paired': True,
                'parse_probs': False,
                'pfix': 'dysgu_reads',
                'pl': 'pe',
                'procs': 1,
                'reference': None,
                'regions': None,
                'regions_mm_only': False,
                'regions_only': False,
                'remap': True,
                'sd': 0.8,
                'search_depth': 20,
                'sites': None,
                'sites_pass_only': True,
                'sites_prob': 0.6,
                'spd': 0.3,
                'sv_aligns': None,
                'symbolic_sv_size': -1,
                'template_size': '',
                'thresholds': {'DEL': 0.45,
                                'INS': 0.45,
                                'INV': 0.45,
                                'DUP': 0.45,
                                'TRA': 0.45},
                'trust_ins_len': True,
                'verbosity': 2,
                'working_directory': 'tempfile',
                'z_breadth': 2,
                'z_depth': 2
            }
    return args


def dysgu_preset_args(mode):
    """
    Returns the default arguments used by dysgu

    :param mode: Choose the 'mode', one of pe | pacbio-sequel2 | nanopore-r9 | nanopore-r10
    :type path: str

    :return: A dict of available arguments
    :rtype: dict
    
    Load a vcf file from dysgu

    :param path: The path to the vcf file
    :type path: str
    :param drop_na_columns: Drop columns that are all NAN
    :type drop_na_columns: bool
    :return: A dataframe of SVs
    :rty
    """
    args = dysgu_default_args()
    presets = {"nanopore-r9": {"mq": 1,
                        "min_support": "auto",
                        "dist_norm": 900,
                        "max_cov": 150,
                        "pl": "nanopore",
                        "remap": "False",
                        "clip_length": -1,
                        "trust_ins_len": "False",
                        "sd": 0.6,
                        "symbolic_sv_size": 50000,
                        "divergence": "auto"
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
                        'thresholds': {'DEL': 0.35,
                                'INS': 0.35,
                                'INV': 0.35,
                                'DUP': 0.35,
                                'TRA': 0.35},
                        "symbolic_sv_size": 50000
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
                      "symbolic_sv_size": 50000
                      },
               "pacbio-revio": {"mq": 1,
                      "min_support": "auto",
                      "dist_norm": 600,
                      "max_cov": 150,
                      "pl": "pacbio",
                      "remap": "False",
                      "clip_length": -1,
                      "trust_ins_len": "True",
                      "sd": 0.4,
                      'thresholds': {'DEL': 0.25,
                                'INS': 0.25,
                                'INV': 0.25,
                                'DUP': 0.25,
                                'TRA': 0.25},
                      "symbolic_sv_size": 50000
                      },
           }
    if mode in presets:
        for k, v in presets[mode].items():
            args[k] = v
    elif mode != "pe":
        raise ValueError("mode must be either: pe | pacbio-sequel2 | nanopore-r9 | nanopore-r10")
    return args



def load_dysgu_vcf(path, drop_na_columns=True):
    """
    Load a vcf file from dysgu

    :param path: The path to the vcf file
    :type path: str
    :param drop_na_columns: Drop columns that are all NAN
    :type drop_na_columns: bool
    :return: A dataframe of SVs
    :rtype: pandas.DataFrame
    """
    df, header, n_fields = vcf_to_df(path)
    if "SVMETHOD" in df:
        del df["SVMETHOD"]
    if drop_na_columns:
        return df.dropna(axis=1, how='all')
    return df


def merge_dysgu_df(*dataframes, merge_distance=500, pick_best=True, add_partners=True):
    """
    Merge calls from dysgu. Input is one or more dataframes with dysgu calls.

    :param dataframes: The input dataframes of dysgu calls to merge
    :type dataframes: pandas.DataFrame
    :param merge_distance: The merging distance, SVs closer than this spacing will be candidates for merging
    :type merge_distance: int
    :param pick_best: A single best SV is chosen for each cluster
    :type pick_best: bool
    :param add_partners: Add information to the output detailing which SVs were merged
    :type add_partners: bool
    :return: The merged data
    :rtype: pandas.DataFrame
    """
    aggressive_ins_merge = True
    table_names = []
    for idx, df in enumerate(dataframes):
        table_names += [idx] * len(df)
    df = pd.concat(dataframes)
    df["table_name"] = table_names
    n_samples = len(dataframes)

    df.reset_index(inplace=True)

    df["event_id"] = df.index
    df["contig"] = df["contigA"]
    df["contig2"] = df["contigB"]

    # preciseA is not saved in vcf, so assume true
    if "preciseA" not in df:
        df["preciseA"] = [1] * len(df)
        df["preciseB"] = [1] * len(df)
    else:
        df["preciseA"].fillna(1, inplace=True)
        df["preciseB"].fillna(1, inplace=True)
    df["contigA"] = [None if pd.isna(i) else i for i in df["contigA"]]
    df["contigB"] = [None if pd.isna(i) else i for i in df["contigB"]]

    potential = [dotdict(set_numeric(i)) for i in df.to_dict("records")]

    bad_i = set([])  # These could not be merged at sample level, SVs probably too close?

    found = merge_events(potential, merge_distance, {}, try_rev=False, pick_best=pick_best, add_partners=add_partners,
                         aggressive_ins_merge=aggressive_ins_merge,
                         same_sample=False)
    ff = defaultdict(set)

    for f in found:
        if f.partners is None:
            ff[f.event_id] = set([])
        else:
            # Remove partners from same sample
            current = f["table_name"]
            targets = set([])
            passed = True

            for item in f["partners"]:  # Only merge with one row per sample
                t_name = df.loc[item]["table_name"]
                if (t_name != current and t_name not in targets and len(targets) < n_samples) or aggressive_ins_merge:
                    targets.add(t_name)
                else:
                    # Merged with self event. Can happen with clusters of SVs with small spacing
                    # e.g. a is merged with b and c, where a is from sample1 and b and c are from sample2
                    # safer not to merge? otherwise variants can be lost
                    passed = False
            if passed:  # enumerate support between components
                g = f["partners"] + [f["event_id"]]
                for t1 in g:
                    for t2 in g:
                        if t2 != t1:
                            ff[t1].add(t2)

    df = df.drop(bad_i)
    if not pick_best:
        if add_partners:
            df["partners"] = [[(df.loc[j].table_name, df.loc[j].event_id) for j in ff[i]] if i in ff else set([]) for i in df.index]
        return DysguSV._mung_df(df)
    df2 = pd.DataFrame.from_records(found)
    if add_partners:
        df2["partners"] = [[(df.loc[j].table_name, df.loc[j].event_id) for j in ff[i]] if i in ff else set([]) for i in df2.index]
    else:
        df2["partners"] = [None] * len(df2)
    return DysguSV._mung_df(df2)



class DysguSV:
    """
    This class is the main interface for calling structural variants using dysgu. To initialize DysguSV,
    provide a reference genome and bam file, via the pysam library. It is also recommended to provide a `sample_name`
    which will be used when saving results to a vcf file.

    :param ref_genome: A reference genome object from the pysam library
    :type ref_genome: pysam.FastaFile
    :param bam: An alignment file from the pysam library
    :type bam: pysam.AlignmentFile
    :param sample_name: The sample name to use in the vcf file
    :type sample_name: str
    :param kwargs: Key-word arguments to modify default options of dysgu
    :type kwargs: dict

    .. highlight:: python
    .. code-block:: python

        import pysam
        from dysgu import DysguSV

        # open a reference genome and alignment file using pysam
        bam = pysam.AlignmentFile('test.bam', 'rb')
        genome = pysam.FastaFile('ref.fasta')

        # Initialise dysgu
        dysgu = DysguSV(genome, bam)

        # Call SVs at a genomic location
        df = dysgu("chr1:1-1000000")

    To options can be provided during initialisation as key-word arguments:

    .. highlight:: python
    .. code-block:: python

        dysgu = DysguSV(genome, bam, min_support=5, mq=20)

    """
    def __init__(self, ref_genome, bam, sample_name="sample", **kwargs):

        if isinstance(ref_genome, str):
            ref_genome = pysam.FastaFile(ref_genome)
        elif not isinstance(ref_genome, pysam.libcfaidx.FastaFile):
            raise ValueError("Expected ref_genome to be of type 'str' or 'pysam.libcfaidx.FastaFile'")
        if not isinstance(bam, pysam.libcalignmentfile.AlignmentFile):
            raise ValueError("Input bam must be a pysam AlignmentFile")
        if bam.is_bam:
            self.kind = "bam"
        elif bam.is_cram:
            self.kind = "cram"
        elif bam.is_sam:
            self.kind = "sam"
        else:
            raise ValueError("Only bam/cram/sam is supported")

        self.ref_genome = ref_genome
        self.bam = bam
        self.sample_name = str(sample_name)
        args = dysgu_default_args()
        if kwargs:
            for k, v in kwargs.items():
                if k not in args:
                    raise ValueError(f"Option {k} not available, see dysgu.default_args().keys() for a full list")
                args[k] = v
        if args["working_directory"] == "tempfile":
            args["working_directory"] = tempfile.mkdtemp()
        self.args = self._fix_args(args)

    def set_option(self, option, value=None):
        """
        Change option(s) for dysgu.

        :param option: The name of the option
        :type option: str
        :param value: The value of the option
        :type value: object
        :return: None
        :rtype: None

        .. highlight:: python
        .. code-block:: python

            dysgu.set_option("min_support", 10)

            # Or provide a mapping of arguments:
            dysgu.set_option({"min_support": 10, "mq": 20, "min_size": 100})

        """
        if isinstance(option, dict):
            for option, value in option.items():
                if option not in self.args:
                    raise ValueError(f"Option {option} not available, see dysgu.default_args().keys() for a full list")
                self.args[option] = value
        elif isinstance(option, str):
            if option not in self.args:
                raise ValueError(f"Option {option} not available, see dysgu.default_args().keys() for a full list")
            self.args[option] = value
        else:
            ValueError("Option must be instance of str or dict")
        self.args = self._fix_args(self.args)

    def _fix_args(self, args):
        # covert some bool's to string format i.e. True to 'True'
        str_format = {'remap', 'drop_gaps', 'trust_ins_len', 'sites_pass_only', 'paired', 'contigs', 'diploid',
                      'regions_only', 'parse_probs', 'all_sites', 'regions_mm_only', 'add_kind', 'verbosity', 'merge_within'}
        args["no_gt"] = args["add_gt"]
        return {k: v if k not in str_format else str(v) for k, v in args.items()}

    def _bed_region_iter(self, regions):
        for item in regions:
            chrom = item[0]
            start = int(item[1])
            end = int(item[2])
            for r in self.bam.fetch(chrom, start, end):
                yield r

    def call_bed_regions(self, regions):
        """
        Call SVs from a list of bed regions. Note,
        bed regions should be sorted by genome starting position, and be non-overlapping. To create a suitable input
        for this function see dysgu.load_bed() and dysgu.merge_intervals() functions.

        :param regions: An iterable of bed regions
        :type regions: iterable
        :return: Dataframe of called SVs
        :rtype: pandas.DataFrame or None

        .. highlight:: python
        .. code-block:: python

            from dysgu import load_bed, merge_intervals

            # load, merge and sort intervals
            bed = load_bed('test.bed')
            bed = merge_intervals(bed, srt=True)

            dysgu = DysguSV(ref, bam)
            df = dysgu.call_bed_regions(bed)

        """
        if not regions:
            return None
        return self.__call__(self._bed_region_iter(regions))

    @staticmethod
    def _mung_df(df):
        df.reset_index(inplace=True)
        df.rename(columns={"contig": "contigA", "contig2": "contigB"}, inplace=True)
        order = ['chrA', 'posA', 'event_id', 'ref_seq', 'variant_seq', 'filter',
                 'sample', 'svtype', 'posB', 'chrB', 'grp_id', 'n_in_grp', 'join_type',
                 'cipos95A', 'cipos95B', 'svlen', 'contigB', 'kind', 'rep', 'rep_sc',
                 'gc', 'n_expansion', 'stride', 'exp_seq', 'ref_poly_bases',
                 'query_overlap', 'su', 'spanning', 'pe', 'supp', 'sc', 'bnd',
                 'svlen_precise', 'type', 'GT', 'GQ', 'NMpri', 'NMsupp', 'NMbase',
                 'MAPQpri', 'MAPQsupp', 'NP', 'maxASsupp', 'sqc', 'scw',
                 'clip_qual_ratio', 'block_edge', 'raw_reads_10kb', 'mcov', 'linked',
                 'neigh', 'neigh10kb', 'ref_bases', 'plus', 'minus', 'strand_binom_t',
                 'n_gaps', 'n_sa', 'n_xa', 'n_unmapped_mates', 'double_clips',
                 'remap_score', 'remap_ed', 'bad_clip_count', 'fcc', 'n_small_tlen',
                 'ras', 'fas', 'inner_cn', 'outer_cn', 'compress', 'ref_rep', 'jitter',
                'contigA', 'right_ins_seq', 'left_ins_seq', 'partners', 'prob']

        df["contigA"].replace(np.nan, "", inplace=True)
        df["contigB"].replace(np.nan, "", inplace=True)

        return df[order]

    def __call__(self, region, sort_df=True):
        """
        Call SVs using dysgu

        :param region: The genomic region to call SVs from
        :type region: str or pysam Iterator
        :param sort_df: Sort the retured dataframe
        :type sort_df: bool
        :return: Dataframe of called SVs
        :rtype: pandas.DataFrame

        .. highlight:: python
        .. code-block:: python

            dysgu = DysguSV(ref, bam)
            df = dysgu("chr1:10000-50000")

            # Using a pysam iterator
            df = dysgu(bam.fetch("chr1", 0, 500000))

        """
        if isinstance(region, str):
            iterable = self.bam.fetch(region=region)
        else:
            iterable = region

        regions = None
        events, site_adder = pipe1(self.args, self.bam, self.kind, regions, self.args['ibam'],
                                                   self.ref_genome, self.sample_name, bam_iter=iterable)
        if len(events) == 0:
            return None
        unused = {"contig2_lc", "contig2_left_weight", "contig2_rc", "contig2_ref_end", "contig2_ref_start",
                  "contig2_right_weight", "contig_lc", "contig_left_weight", "contig_rc", "contig_ref_end",
                  "contig_ref_start", "contig_right_weight", "modified", "preciseA",
                  "preciseB", "query_gap", "remapped", "site_info"}  # "partners",
        df = pd.DataFrame.from_records([{k: v for k, v in to_dict(e).items() if k not in unused} for e in events])

        df = self.apply_model(df)
        args = self.args
        if args["sites"]:
            df = post_call_metrics.update_prob_at_sites(df, events, args["thresholds"],
                                                        parse_probs=args["parse_probs"] == "True",
                                                        default_prob=args["sites_prob"])
            df["site_id"] = ["." if not s else s.id for s in df["site_info"]]
            if args["all_sites"] == "True":
                raise NotImplementedError("all-sites is not supported using the python-api currently")

        if sort_df:
            df = df.sort_values(["chrA", "posA", "event_id"])

        df["sample"] = [self.sample_name] * len(df)

        # fix variant seq column
        df["variant_seq"] = [i.upper() if svt == "INS" and i is not None and len(i) > 0 else f"<{svt}>" for i, svt in zip(df["variant_seq"], df["svtype"], strict=False)]
        return self._mung_df(df)

    def to_vcf(self, dataframe, output_file):
        """
        Save dysgu SV calls to a vcf file

        :param dataframe: A dataframe of called SVs from dysgu
        :type dataframe: pandas.DataFrame
        :param output_file: The file handle to write the vcf file to
        :type output_file: file
        :return: None
        :rtype: None

        .. highlight:: python
        .. code-block:: python

            with open(path, "w") as out:
                dysgu.to_vcf(passed, out)

        """
        contig_header_lines = ""
        for item in self.bam.header["SQ"]:
            contig_header_lines += f"\n##contig=<ID={item['SN']},length={item['LN']}>"
        args = self.args
        to_vcf(dataframe, args, {self.sample_name}, output_file, show_names=False, contig_names=contig_header_lines)

    def apply_model(self, df):
        """
        Apply a machine leaning model to the dataframe. The model configuration is determined by the options set on
        the DysguSV class. For example, to use a non-diploid model, first set `diploid=False`:

        :param df: The input dataframe to apply the machine learning model to
        :type df: pandas.DataFrame
        :return: Dataframe with a modified 'prob' column
        :rtype: pandas.DataFrame

        .. highlight:: python
        .. code-block:: python

            dysgu = DysguSV(ref, bam)
            dysgu.set_option("diploid", False)
            df = dysgu("chr1:10000-50000")

        """
        args = self.args
        df = post_call_metrics.apply_model(df, args["pl"], args["contigs"], args["diploid"], args["thresholds"])
        return df

