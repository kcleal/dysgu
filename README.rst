
.. image:: dysgu/logo.png
    :align: left

.. |Generic badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
   :target: http://bioconda.github.io/recipes/dysgu/README.html

.. |Li badge| image:: https://anaconda.org/bioconda/dysgu/badges/license.svg
   :target: https://github.com/kcleal/dysgu/blob/master/LICENSE.md

dysgu (pronounced *duss-key*) is a set of command line tools and `python-API <https://kcleal.github.io/dysgu/API.html>`_,
for calling structural variants using paired-end or long read sequencing data. See recent long-read benchmarks `here <https://github.com/kcleal/SV_Benchmark_CMRG>`_, and `here <https://github.com/kcleal/SV_benchmark_PacBio_HiFi>`_.

|Generic badge| |Li badge|

`⚙️ Installation`_

`🚀 Quick start`_

`🎯 Calling SVs`_

`🚦Filtering SVs`_

`➕ Merging SVs`_

`♋ Somatic SVs / tumor-normal calling / pool-of-normals`_

`🔍 Genotype list of sites`_

`🔪 Regions of interest / excluding regions`_

`🔧 Useful parameters`_

`🚑 Issues`_

`🐍 Python API`_

`🎓 Citation`_

----

⚙️ Installation
---------------

Dysgu can be installed via pip using Python >=3.10 and has been tested on linux and MacOS.
The list of python packages needed can be found in requirements.txt.
To install::

    pip install dysgu

Or, from conda::

    conda install -c conda-forge -c bioconda dysgu


To build from source, run the install script: ``bash INSTALL.sh``.

Alternatively, pull from `dockerhub <https://hub.docker.com/repository/docker/kcleal/dysgu/>`_::

    docker pull kcleal/dysgu

Run tests::

    $ dysgu test

🚀 Quick start
--------------
Available commands::

    dysgu run              # Run using default arguments, wraps fetch and call commands
    dysgu fetch            # Separate SV reads from input bam file
    dysgu call             # SV calling
    dysgu merge            # Merge calls from multiple samples
    dysgu filter           # Filter SVs, find somatic SVs (version >= 1.5.0)
    dysgu test             # Run basic tests

For help use::

    dysgu --help
    dysgu command --help

To use the python-API see the `documentation <https://kcleal.github.io/dysgu/API.html>`_,
or `jupyter notebook <https://github.com/kcleal/dysgu/blob/master/dysgu_api_demo.ipynb>`_,


🎯 Calling SVs
--------------

Paired-end reads
****************
Dysgu was developed to work with paired-end reads aligned using bwa mem, although other aligners will work.
To call SVs, a sorted and indexed .bam/cram is needed plus an indexed reference genome (fasta format).
Also a working directory must be provided to store temporary files.
There are a few ways to run dysgu depending on the type of data you have.
For paired-end data the `run` command is recommended which wraps `fetch` and `call`::

    dysgu run reference.fa temp_dir input.bam > svs.vcf

This will first run the `fetch` command that will create a temporary bam file plus other analysis files
in the directory `temp_dir`. These temporary files are then analysed using the `call` program.
To make use of multiprocessing, set the "-p" parameter::

    dysgu run -p4 reference.fa temp_dir input.bam > svs.vcf

Dysgu also accepts reads from stdin. In this example, the --clean flag will remove temp files on completion::

    samtools view -bh samp.bam chr1:0-1000000 | dysgu run --clean reference.fa temp_dir - > svs.vcf

Long reads
**********
Dysgy is designed to work with long reads aligned using minimap2 or ngmlr. Use the 'call' pipeline if
starting with a bam file, or 'run' if starting with a cram::

    dysgu call --mode pacbio-revio reference.fa temp_dir input.bam > svs.vcf

    dysgu call --mode nanopore-r10 reference.fa temp_dir input.bam > svs.vcf

Presets are available using the `--mode` option for PacBio `pacbio-sequel2 | pacbio-revio`,
and ONT `nanopore-r9 | nanopore-r10`.

Since v1.8, higher calling accuracy can be achieved by adding "haplotags" to your input bam/cram file
(using hiphase / whatshap / longshot etc). Dysgu reads HP and PS tags automatically and produces
phased output calls.

If you are using using reads with higher error rates, or are unsure of the read-accuracy,
it is recommended to set '--divergence auto', to infer a more conservative sequence divergence.
The default is set at 0.02 which can be too stringent for lower accuracy reads and will result in
more reads being filtered and lower sensitivity::

    dysgu call --divergence auto --mode nanopore reference.fa temp_dir input.bam > svs.vcf

Pipeline overview
~~~~~~~~~~~~~~~~~
The first stage of the "run" pipeline is to separate SV-associated reads - split/discordant reads,
and reads with a soft-clip >= clip_length (15 bp by default for paired-end reads).
This is achieved using the `fetch` command which can be run independently if needs be::

    dysgu fetch samp1_temp input.bam

All SV associated reads will be placed in `samp1_temp/input.dysgu_reads.bam`.
The next stage of the pipeline is to call SVs using the `call` command. Additionally, the `--ibam` option is recommended for paired-end data so dysgu can infer insert
size metrics from the main alignment file. If this is not provided, dysgu will use the input.bam in the samp1_temp folder which may be less accurate. Alternatively,
the insert size can be specified manually using the -I option::

    dysgu call --ibam all_reads.bam reference.fa temp_dir temp_dir/temp_dir.dysgu_reads.bam > svs.vcf

Models available
~~~~~~~~~~~~~~~~~
There are a choice of three models per read type. By default, a diploid model will be used that takes into account
changes in read-depth around break sites. This model is
preferred as it often attains higher precision in germline whole-genome samples. However, for somatic samples (e.g. tumors) copy
number changes, poly-clonality or poly-ploidy can lead to events with low allelic fraction. For such samples, a non-diploid
model might work better. This is selected by applying `--diploid False`. A model with no information on allelic fraction
will then be utilized.

Finally, if the diploid/non-diploid models are not picking up your SV of interest, a simpler model can be used with the
`--contigs False` option. This model has all sequence-related metrics removed, so only read-support information is
retained. In general the performance of models follows diploid > non-diploid > no-contigs.

It is also possible to switch models post-calling using the python-API. For an example of how to do this,
see the dysgu_api_demon.ipynb

Resource requirements
~~~~~~~~~~~~~~~~~~~~~
Using a single core and depending on hard-drive speed, dysgu usually takes ~1h to analyse a 30X coverage genome of 150 bp paired-end reads and
uses < 6 GB memory. Also note that when `fetch` is utilized (or using run command), a large temp file is generated consisting of SV-associated reads >5 Gb in size.


🚦Filtering SVs
----------------
The filtering command is quite flexible and can be used to filter a single sample,
or filter against a normal sample, or a panel of normals/cohort.
If the filter command is used with a single input vcf (no normals or cohort), filtering will remove lower quality events::

    dysgu filter input.vcf > output.vcf

Filtering is generally recommended after any merging has been performed, or if you are analysing only a single sample.

If a normal vcf is supplied, then input calls will be removed if they overlap with events in the normal vcf::

    dysgu filter --normal-vcf normal.vcf input.vcf > output.vcf

Additionally, you can provide bam files to filter against. This will make the filtering much more stringent as each
alignment file you provide will be checked for reads that match your input calls. If supporting reads are found then
the input call will be removed. Note, this also makes filtering much slower. For large cohorts a random sample of
bams can be used for filtering using the `--random-bam-sample Int` option::

    dysgu filter input.vcf normal.bam > output.vcf  # normal bam only
    dysgu filter --normal-vcf normal.vcf  input.vcf normal.bam > output.vcf

Dysgu will understand the sample-name in vcf and bam files, so if you use "*.bam" syntax, then the input sample will
not be used for filtering.

Other filtering option are detailed below.

Remove events with low probability::

    dysgu filter --min-prob 0.2 input.vcf > output.vcf

Remove events with low support fraction::

    dysgu filter --support-fraction 0.15 input.vcf > output.vcf

Re-label events with probability >= 0.3 as PASS::

    dysgu filter --pass-prob 0.3 input.vcf > output.vcf

Use normal bams to filter common/germline structural variants::

    dysgu filter input.vcf normals/*.bam > output.vcf
    dysgu filter input.vcf list_of_normals.txt > output.vcf


➕ Merging SVs
--------------
If you plan on merging samples, it is recommended that the '-v2' option be used when running the 'run/call' modules; this will
ensure that all consensus sequences will be reported in the vcf file to help with downstream merging.
Multiple output vcf files can be merged, e.g. tumor.vcf and normal.vcf, or illumina.vcf and pacbio.vcf::

    dysgu merge *.vcf > combined.vcf

For large numbers of samples, an input list can be used, and merging can be performed in parallel (by chromosome and SV type)::

    dysgu merge -p24 --input-list samples.txt --wd wd > combined.vcf

Merging SVs between platforms at multiallelic/complex sites is still tricky and there is a trade off between under merging
(leading to duplication) and over merging (leading to loss of multiallelic/complex SVs). Setting the '--merge-within True' option will perform
a single round of merging for each input file before merging across input files. This will shift the balance to over merging, but reduces the
problem of duplication::

    dysgu merge --merge-within True pacbio.vcf illumina.vcf > combined.vcf


♋ Somatic SVs / tumor-normal calling / pool-of-normals
------------------------------------------------------

For tumor/normal pairs, the recommended workflow is to call SVs independently in each sample, then obtain tumor specific (somatic) SVs by running dysgu filter::

    dysgu run ref.fa wd_t tumour.bam > tumor.vcf
    dysgu run ref.fa wd_n normal.bam > normal.vcf
    dysgu filter --normal-vcf normal.vcf tumour.vcf normal.bam > somatic.vcf

The output vcf will contain SVs that are deemed to be unique in the tumor sample.

Unique SV can also be identified when compared to a cohort vcf or list of bam files. A third-party vcf of common SVs can be used (provided 'SVTYPE' is listed in the info column). Or,
cohort SVs can be merged using `dysgu merge`, before filtering to get unique SVs::

    dysgu merge *.vcf > merged.vcf
    dysgu filter --normal-vcf merged.vcf sample1.vcf *.bam > sample1_unique.vcf
    dysgu filter --normal-vcf merged.vcf sample1.vcf list_of_normals.txt > sample1_unique.vcf

Here, sample1.vcf and merged.vcf can contain multiple samples, although if sample1.vcf is multi-sample, you must provide '--target-sample' to indicate which sample to filter.
The output sample1_somatic.vcf will be a single sample vcf containing unique SVs.

Sample names are respected from the vcf and bam file headers (or filenames), so `sample1` will be ignored from the normal-vcf and list of bams.
To keep all SVs in the output, use ``--keep-all``. Filtered SVs will be labelled 'normal', 'lowProb' or 'lowSupport' in the filter column.

Increasing the number of bams to filter against will slow down filtering, but should increase specificity. To set a
limit on the number of bams to filter against, a random sample can be drawn from the input list,
e.g. draw 5 random bam samples from the input list to filter against using::

    dysgu filter --random-bam-sample 5 --normal-vcf merged.vcf sample1.vcf *.bam


Also a target VCF can be filtered against a normal vcf if desired (without alignment files)::

    dysgu filter --normal-vcf normal.vcf sample1.vcf

By default, SV calls with a PROB value < ``--min-prob`` are removed from the final output,
and SV calls with a PROB value >= ``--pass-prob`` will be re-labelled as PASS in the output. However, these
thresholds currently require tuning depending on sequencing platform, coverage and the size of the cohort used for filtering.
Suitable values for `--pass-prob` often lie in the range 0.2 - 0.4. For paired-end reads, a pass-prob of around 0.35 can work well, whereas for long-reads a lower threshold of 0.2 can work better e.g::

    dysgu filter --pass-prob 0.2 --min-prob 0.1 --normal-vcf normal.vcf tumour.vcf normal.bam > somatic.vcf

To quickly test and visualise different filtering thresholds, output can be piped to the command line tool `GW <https://github.com/kcleal/gw>`_, which will display the results to screen for inspection::

    dysgu filter --pass-prob 0.2 filtered.vcf | \
    gw hg38 -b normal.bam -b tumor.bam -v -


🔍 Genotype list of sites
-------------------------
Calls from multiple samples can be merged into a unified site list::

    dysgu run -v2 ref.fa wd1 sample1.bam > sample1.vcf
    dysgu run -v2 ref.fa wd2 sample2.bam > sample2.vcf
    dysgu merge sample1.vcf sample2.vcf > merged.vcf

This list can be used to re-genotype at the sample level. Here, to save time, the temporary files in the working directory 'wd1' are re-used::

    dysgu call --ibam sample1.bam --sites merged.vcf ref.fa wd1 wd1/sample1.dysgu_reads.bam > sample1.re_geno.vcf

This is equivalent to running::

    dysgu run --sites merged.vcf ref.fa wd1 sample1.bam > sample1.re_geno.vcf

Dysgu can also accept --sites from other sources, for example calls from other SV callers or read-types can be provided::

    dysgu run --sites manta.diploidSVs.vcf ref.fa wd sample1.bam > sample1.vcf

This can help discovery of events with low read-support.

To output all variants in --sites including those with genotype 0/0 in the input sample, set '--all-sites True'.

By default if a matching call is found in both --sites and the input sample, then the probability value
(PROB value in the FORMAT field of the output vcf) of the call will be modified. This behavior can be controlled by setting the
--sites-prob option (default value is 0.6), controlling the probability that a matching call in --sites is a true
variant in the input sample. To turn this behavior off, set the --sites-prob value to 0.5, which implies an even chance that a matching site
in --sites is also a true variant in the input sample. For related individuals or samples, or if the
--sites are from a trusted source, a higher --sites-prob value is recommended e.g. --sites-prob 0.8.

If the --sites vcf file is from a previous dysgu run, the PROB values can be utilized by setting '--parse-probs True'. This
option can work well when using dysgu calls from a related individual.

Also of note, the ``--ignore-sample-sites`` option is set to True by default. This results in the input sample name (from the bam SM tag)
 being ignored from a multi-sample sites file. This may not be the deired behavior if trying to re-genotype a sample using different
 read types, for example.


🔪 Regions of interest / excluding regions
------------------------------------------
Regions of the genome can be skipped from analysis by providing a .bed file using the `--exclude` option. This option
takes precedence over the options detailed below, and acts as a hard filter, removing regions of the genome from analysis.

Dysgu provides two ways to analyse regions of interest. Target genomic regions can be specified using a .bed file with
the --search option. This will also act as a hard filter, limiting analysis only to those regions, while regions outside
will be ignored.

Alternatively, regions can be specified using the --regions option (.bed file). If this option is used, all reads not
excluded by the --exclude/--search options will be analysed. Variants will then be
labelled in the output vcf according to their intersection with those regions. The INFO > KIND column will be labelled
with either 'intra-regional' - both SV ends within same interval, 'extra-regional' - neither SV end in an interval,
'inter-regional' - SV ends in separate intervals, or 'hemi-regional' - one SV end in an interval. These labels may be
useful for some targeted sequencing experiments.

Additionally, there is also the --regions-only option. The option is only available for 'dysgu call'. If this is set to 'True', then dysgu will search all reads in
--regions and also analyse any mate-pairs that do not overlap those regions of interest. This method can be quicker to
run when the regions of interest are small relative to the genome. However, this option can consume a lot of memory if the
regions are large, so use with caution.

For deep targeted sequencing experiments, the --regions-mm-only option can also be used, which can help prevent over
clustering of reads. When set to 'True', dysgu will only use minimizer based clustering within the intervals specified
by --regions.

Also of note, it is possible to use --exclude, --search, and --regions at the same time.


🔧 Useful parameters
--------------------
The most important parameter affecting sensitivity is --min-support, lower values increase sensitivity but also runtime.

The --max-cov parameter may need to be adjusted for high coverage samples (default is 200), or samples that might have
high copy number aberrations. Only reads with mapq >= `--mq` threshold count towards coverage values and regions with coverage exceeding `max-cov` are ignored for SV calling.
Dysgu can automatically infer a max-cov value for bam files by setting `--max-cov auto`, which
will correspond to ~6*whole-genome-coverage by default. However using 'auto', is only recommended for whole-genome samples.
A helper script can be used to suggest different max-cov values with respect to mean genome coverage, for example
to use of threshold of 25 x mean genome coverage::


    max_cov=$(python suggest_max_coverage.py -y 25 input.bam)
    >>> Read-length 148.0 bp, mean whole-genome coverage estimate: 31.88, max-cov ~ 797

    dysgu run --max-cov $max_cov reference.fa temp_dir input.bam > svs.vcf

The --thresholds parameter controls the probability value at which events are labelled with a
'PASS', increasing these values increases precision at the expense of sensitivity.

The verbosity of contig reporting can be controlled using '-v/--verbosity'. If you plan to use "merge" on output files,
it is a good idea to use "-v2" as contig sequences can help with merging.

--trust-ins-len applies to long-read data (pacbio, nanopore modes). If set to 'True', insertion length as stated in
the alignment cigar string is assumed to be correct and more stringent clustering is utilized. This can improve sensitivity at multi-allelic
sites but at the expense of increasing duplicate true-positive calls that arise mostly at SVs with
ambiguous candidate alignments.

--divergence applies to long reads only, and measures the proportion of non-reference cigar operations (deletions, insertions)
compared to matching reference bases. Reads that have anomalous divergence at the ends of the read are ignored during calling.


🚑 Issues
---------
- Currently cram files are only supported when using the "run" command. This is because pysam cannot use seek on a cram file.

- If the temp file created during the fetch stage of the pipeline is too big, the --compression level can be set to reduce space.

- If dysgu is taking a long time to run, this could be due to the complexity of the sample. Dysgu will try and generate contigs from clusters of soft-clipped reads and remap these to the reference genome. In this case consider increasing the `clip-length` or setting `--contigs False`, or `--remap False`. Alternatively you might need to check your sample for anomalous sequences and adapter content.

- If dysgu is consuming a large amount of memory, you can try the --low-mem flag.

- If sensitivity is lower than expected, check that the insert size was inferred accurately (provide manually using the `-I`), and divergence is set appropriately.

- If you input data or aligner do not seem to be working well with dysgu, please get in touch clealk@cardiff.ac.uk


🐍 Python API
-------------

Dysgu can also be used from a python script. A full demo of the API can be found in the
`ipython notebook <https://github.com/kcleal/dysgu/blob/master/dysgu_api_demo.ipynb>`_,. In this example, dysgu is
used to call SVs on the first 10 Mb of chr1:

.. code-block:: python

    import pysam
    from dysgu import DysguSV

    # open input bam and reference file
    bam = pysam.AlignmentFile('sample.bam', 'rb')
    genome = pysam.FastaFile('ucsc.hg19.fasta')

    # initiate dysgu
    dysgu = DysguSV(genome, bam)

    # call SVs (results will be a pandas dataframe)
    results = dysgu(bam.fetch('chr1', 0, 10_000_000))

    # after analysis, save to a vcf file
    with open("output.vcf", "w") as out:
        dysgu.to_vcf(results, out)

The API can also be used to apply different machine-learning models, merge SVs, and call SVs using target bed regions.

🎓 Citation
-----------
To cite dysgu, or to learn more about implementation details please see:

https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac039/6517943



