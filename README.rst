
========
dysgu-SV
========

dysgu (pronounced *duss-key*) is a set of command line tools and python-API,
for calling structural variants using paired-end or long read sequencing data.


Installation
------------
Dysgu requires Python >=3.7 - 3.10 and has been tested on linux and MacOS.
The list of python packages needed can be found in requirements.txt.
To install::

    pip install numpy dysgu

To build from source, a >=c++11 compiler is needed. For convenience use the install script::

    pip install numpy
    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu
    bash INSTALL.sh

Or manually::

    pip install numpy
    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu/dysgu/htslib
    make
    cd ../../
    pip install -r requirements.txt
    pip install .

Pull from `dockerhub <https://hub.docker.com/repository/docker/kcleal/dysgu/>`_::

    docker pull kcleal/dysgu

Run tests::

    $ dysgu test


Usage
-----
Available commands::

    dysgu run     # Run using default arguments, wraps fetch and call commands
    dysgu fetch   # Separate SV reads from input bam file
    dysgu call    # SV calling
    dysgu merge   # Merge calls from multiple samples
    dysgu test    # Run basic tests

For help use::

    dysgu --help
    dysgu command --help


Calling SVs
~~~~~~~~~~~

Paired-end reads
****************
Dysgu was developed to work with paired-end reads aligned using bwa mem, although other aligners will work. To call SVs, a sorted and indexed .bam/cram is needed plus an indexed reference genome (fasta format). Also a working directory must
be provided to store temporary files. There are a few ways to run dysgu depending on the type of data you have.
For paired-end data the `run` command is recommended which wraps `fetch` and `call`::

    dysgu run reference.fa temp_dir input.bam > svs.vcf

This will first run the `fetch` command that will create a temporary bam file plus other analysis files in the directory `temp_dir`. These temporary files are then analysed using the `call` program.
To make use of multiprocessing, set the "-p" parameter::

    dysgu run -p4 reference.fa temp_dir input.bam > svs.vcf

Dysgu also accepts reads from stdin. In this example, the --clean flag will remove temp files on completion::

    samtools view -bh samp.bam chr1:0-1000000 | dysgu run --clean reference.fa temp_dir - > svs.vcf

Long reads
**********
Dysgy has been designed to work with long reads aligned using minimap2 or ngmlr. For very long reads (Oxford nanopore), the `fetch` stage of the pipeline is not necessary, so the `call` command should be used directly.
For PacBio Sequel II HiFi reads, the `run` command is generally recommended as it results in lower run times although at the expense of generating additional temp files in the working directory::

    dysgu call --mode pacbio reference.fa temp_dir input.bam > svs.vcf

    dysgu call --mode nanopore reference.fa temp_dir input.bam > svs.vcf



Pipeline overview
~~~~~~~~~~~~~~~~~
The first stage of the "run" pipeline is to separate SV-associated reads - split/discordant reads,
and reads with a soft-clip >= clip_length (15 bp by default).
This is achieved using the `fetch` command, which can be run independently if needs be::

    dysgu fetch samp1_temp input.bam


All SV associated reads will be placed in `samp1_temp/input.dysgu_reads.bam`.
The next stage of the pipeline is to call SVs using the `call` command. Additionally, the `--ibam` option is recommended for paired-end data so dysgu can infer insert
size metrics from the main alignment file. If this is not provided, dysgu will use the input.bam in the samp1_temp folder which may be less accurate. Alternatively,
the insert size can be specified manually using the -I option::

    dysgu call --ibam all_reads.bam reference.fa temp_dir temp_dir/temp_dir.dysgu_reads.bam > svs.vcf


Merging SVs from multiple files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you plan on merging samples, it is recommended that the -'v2' option be used when running the 'run/call' modules; this will
ensure that all consensus sequences will be reported in the vcf file to help with downstream merging.
Multiple output vcf files can be merged, e.g. tumor.vcf and normal.vcf, or illumina.vcf and pacbio.vcf::

    dysgu merge sample1.vcf sample2.vcf > combined.vcf

Merging SVs between platforms at multiallelic/complex sites is still tricky and there is a trade off between under merging
(leading to duplication) and over merging (leading to loss of multiallelic/complex SVs). Setting the '--merge-within True' option will perform
a single round of merging for each input file before merging across input files. This will shift the balance to over merging, but reduces the
problem of duplication::

    dysgu merge --merge-within True pacbio.vcf illumina.vcf > combined.vcf


Models available
----------------
There are a choice of three models per read type. By default, a diploid model will be used that takes into account
changes in read-depth around break sites. This model is
preferred as it often attains higher precision in germline whole-genome samples. However, for somatic samples (e.g. tumors) copy
number changes, poly-clonality or poly-ploidy can lead to events with low allelic fraction. For such samples, a non-diploid
model might work better. This is selected by applying `--diploid False`. A model with no information on allelic fraction
will then be utilized.

Finally, if the diploid/non-diploid models are not picking up your SV of interest, a simpler model can be used with the
`--contigs False` option. This model has all sequence-related metrics removed, so only read-support information is
retained. In general the performance of models follows diploid > non-diploid > no-contigs.

Specifying regions of interest / excluding regions
--------------------------------------------------

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

Additionally, there is also the --regions-only option. If this is set to 'True', then dysgu will search all reads in
--regions and also analyse any mate-pairs that do not overlap those regions of interest. This method can be quicker to
run when the regions of interest are small relative to the genome.

For deep targeted sequencing experiments, the --regions-mm-only option can also be used, which can help prevent over
clustering of reads. When set to 'True', dysgu will only use minimizer based clustering within the intervals specified
by --regions.

Also of note, it is possible to use --exclude, --search, and --regions at the same time.


Genotype list of sites
----------------------
Calls from multiple samples can be merged into a unified site list::

    dysgu run -v2 ref.fa wd1 sample1.bam > sample1.vcf
    dysgu run -v2 ref.fa wd2 sample2.bam > sample2.vcf
    dysgu merge sample1.vcf sample2.vcf ... > merged.vcf

This list can be used to re-genotype at the sample level. Here, to save time, the temporary files in the working directory 'wd1' are re-used::

    dysgu call --ibam sample1.bam --sites merged.vcf ref.fa wd1 wd1/sample1.dysgu_reads.bam > sample1.re_geno.vcf

This is equivalent to running::

    dysgu run --sites merged.vcf ref.fa wd1 sample1.bam > sample1.re_geno.vcf

Dysgu can also accept --sites from other sources, for example calls from other SV callers or read-types can be provided::

    dysgu run --sites manta.diploidSVs.vcf ref.fa wd sample1.bam > sample1.vcf

This can especially help discovery of events with low read-support.

To output all variants in --sites including those with genotype 0/0 in the input sample, set '--all-sites True'.

By default if a matching call is found in both --sites and the input sample, then the probability value
(PROB value in the FORMAT field of the output vcf) of the call will be modified. This behavior can be controlled by setting the
--sites-prob option (default value is 0.6), controlling the probability that a matching call in --sites is a true
variant in the input sample. To turn this behavior off, set the --sites-prob value to 0.5, which implies an even chance that a matching site
in --sites is also a true variant in the input sample. For related individuals or samples, or if the
--sites are from a trusted source, a higher --sites-prob value is recommended e.g. --sites-prob 0.8.

If the --sites vcf file is from a previous dysgu run, the PROB values can be utilized by setting '--parse-probs True'. This
option can work well when using dysgu calls from a related individual.


Useful parameters
-----------------
The most important parameter affecting sensitivity is --min-support, lower values increase sensitivity but also runtime.

The --max-cov parameter may need to be adjusted for high coverage samples (default is 200), or samples that might have
high copy number aberrations. Regions with coverage exceeding `max-cov` are ignored for SV calling.
Dysgu can automatically infer a max-cov value for bam files by setting `--max-cov auto`, which
will correspond to ~6*whole-genome-coverage by default. However using 'auto', is only recommended for whole-genome samples.
A helper script can be used to suggest different max-cov values with respect to mean genome coverage, for example
to use of threshold of 25 x mean genome coverage::


    max_cov=$(python scripts/suggest_max_coverage.py -y 25 input.bam)
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


Resource requirements
---------------------
Using a single core and depending on hard-drive speed, dysgu usually takes ~1h to analyse a 30X coverage genome of 150 bp paired-end reads and
uses < 6 GB memory. Also note that when `fetch` is utilized (or using run command), a large temp file is generated consisting of SV-associated reads >5 Gb in size.


Issues
------
- Currently cram files are only supported when using the "run" command. This is because pysam cannot use seek on a cram file.

- If the temp file created during the fetch stage of the pipeline is too big, the --compression level can be set to reduce space.

- If dysgu is taking a long time to run, this could be due to the complexity of the sample. Dysgu will try and generate contigs from clusters of soft-clipped reads and remap these to the reference genome. In this case consider increasing the `clip-length` or setting `--contigs False`, or `--remap False`. Alternatively you might need to check your sample for anomalous sequences and adapter content.

- If dysgu is consuming a large amount of memory, you can try the --low-mem flag.

- If sensitivity is lower than expected for paired-end data, check that the insert size was inferred accurately, and provide manually using the `-I` option otherwise.

- If you input data or aligner do not seem to be working well with dysgu, please get in touch clealk@cardiff.ac.uk


Citation
--------
To cite dysgu, or to learn more about implementation details please see:

https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac039/6517943



