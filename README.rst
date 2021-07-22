
========
dysgu-SV
========

dysgu (pronounced *duss-key*) is a set of tools for calling structural variants using paired-end or long read sequencing data.


Installation
------------
Dysgu requires Python>=3.7 and has been tested on linux and MacOS.
The list of python packages needed can be found in requirements.txt.
To install::

    pip install dysgu

To build from source, a >=c++11 compiler is needed. For convenience use the install script::

    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu
    bash INSTALL.sh

Or manually::

    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu/dysgu/htslib
    make
    cd ../../
    pip install -r requirements.txt
    pip install .

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


For example, to save time, `dysgu fetch` can be run in a stream during mapping/sorting. In the example below, dysgu reads from stdin and
all SV associated reads will be placed in `temp_dir/temp_dir.dysgu_reads.bam`, and all input alignments will be placed in all_reads.bam::

    dysgu fetch -r all_reads.bam temp_dir -

The next stage of the pipeline is to call SVs using the `call` command. Additionally, the `--ibam` option is recommended for paired-end data so dysgu can infer insert
size metrics from the main alignment file. If this is not provided, dysgu will use the input.bam in the samp1_temp folder which may be less accurate. Alternatively,
the insert size can be specified manually using the -I option::

    dysgu call --ibam all_reads.bam reference.fa temp_dir temp_dir/temp_dir.dysgu_reads.bam > svs.vcf


Merging SVs from multiple files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiple output vcf files can be merged, e.g. tumor.vcf and normal.vcf, or illumina.vcf and pacbio.vcf::

    dysgu merge pacbio.vcf illumina.vcf > combined.vcf

For help use::

    dysgu --help
    dysgu command --help


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
it is usually a good idea not to use "-v0" as contig sequences can help with merging.

Resource requirements
---------------------
Using a single core and depending on hard-drive speed, dysgu usually takes ~1h to analyse a 30X coverage genome of 150 bp paired-end reads and
uses < 6 GB memory. Also note that when `fetch` is utilized (or using run command), a large temp file is generated consisting of SV-associated reads >5 Gb in size.


Issues
------
Currently cram files are only supported when using the "run" command. This is because pysam cannot use seek on
a cram file.

If the temp file created during the fetch stage of the pipeline is too big, the --compression level can be
set to reduce space.

If dysgu is taking a long time to run, this could be due to the complexity of the sample.
Dysgu will try and generate contigs from clusters of soft-clipped reads and remap these to the reference genome.
In this case consider increasing the `clip-length` or setting `--contigs False`, or `--remap False`.
Alternatively you might need to check your sample for anomalous sequences and adapter content.

If sensitivity is lower than expected for paired-end data, check that the insert size was inferred accurately, and
provide manually using the `-I` option otherwise.

If you input data or aligner do not seem to be working well with dysgu, please get in touch clealk@cardiff.ac.uk


Citation
--------
To cite dysgu, or to learn more about implementation details please see our pre-print:

https://www.biorxiv.org/content/10.1101/2021.05.28.446147v1.full


