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
    dysgu fetch   # Seperate SV reads from input bam file
    dysgu call    # SV calling
    dysgu merge   # Merge calls from multiple dysgu call runs
    dysgu test    # Run basic tests

Calling SVs
~~~~~~~~~~~

Paired-end reads
****************
Dysgu was developed to work with paired-end reads aligned using bwa mem. To call SVs, a sorted and indexed .bam/cram is needed plus an indexed reference genome (fasta format). Also a working directory must
be provided to store temporary files. There are a few ways to run dysgu depending on the type of data you have.
For paired-end data the `run` command is recommended which wraps `fetch` and `call`::

    dysgu run reference.fa temp_dir input.bam > svs.vcf

This will first call `fetch` which will create a temporary bam file plus other analysis files in the working directory `samp1_temp`. These temporary files are then analysed using the `call` program.

Long reads
**********
Dysgy has been designed with long reads aligned using minimap2 or ngmlr. For very long reads (Oxford nanopore), the `fetch` stage of the pipeline is not necessary, so the `call` command should be used directly.
For PacBio Sequel II HiFi reads, the `run` command is recommended as it results in lower run times although at the expense of generating additional temp files in the working directory::

    dysgu call --mode pacbio reference.fa temp_dir input.bam > svs.vcf
    dysgu call --mode nanopore reference.fa temp_dir input.bam > svs.vcf


Fetching SV reads
~~~~~~~~~~~~~~~~~
To separate SV-associated reads, use `fetch` on an existing .bam file::

    dysgu fetch samp1_temp input.bam


To save time, `dysgu fetch` can be run in a stream during mapping/sorting. Here, dysgu reads from stdin and
all SV associated reads will be placed in `temp_dir/temp_dir.dysgu_reads.bam`, and all input alignments will be placed in all_reads.bam::

    dysgu fetch -r all_reads.bam temp_dir -

SVs can be subsequently called using the `call` command. Additionally, the `--ibam` option is recommended for paired-end data so dysgu can infer insert
size metrics from the main alignment file. If this is not provided, dysgu will use the input.bam in the samp1_temp folder which may be less accurate::

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
preferred as it often attains higher precision in germline samples. However, for somatic samples e.g. tumors, copy
number changes, poly-clonality or poly-ploidy can lead to events with low allelic fraction. For such samples, a non-diploid
model might work better. This is selected by applying `--diploid False`. A model with no information on allelic fraction
will then be utilized.

Finally, if the diploid/non-diploid models are not picking up your SV of interest, a simpler model can be used with the
`--contigs False` option. This model has all sequence-related metrics removed, so only read-support information is
retained. In general the performance of models follows diploid > non-diploid > no-contigs.

Resource requirements
---------------------
Using a single core and depending on hard-drive speed, dysgu usually takes ~1h to analyse a 30X coverage genome of 150 bp paired-end reads and
uses < 6 GB memory. Also note that when `fetch` is utilized (or using run command), a large temp file is generated consisting of SV-associated reads >5 Gb in size.

Issues
------
If dysgu is taking a long time to run, this could be due to the complexity of the sample.
Dysgu will try and generate contigs from clusters of soft-clipped reads and remap these to the reference genome.
In this case consider increasing the `clip-length` or setting `--contigs False`, or `--remap False`.
Alternatively you might need to check your sample for anomalous sequences and adapter content.

If sensitivity is lower than expected for paired-end data, check that the insert size was inferred accurately, and
provide manually using the `-I` option otherwise.

If you input data or aligner do not seem to be working well with dysgu, please get in touch clealk@cardiff.ac.uk
