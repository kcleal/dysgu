========
dysgu-SV
========

[Build Status](https://travis-ci.com/kcleal/dysgu.svg?token=ggp1k8nRaRrwARctVfix&branch=master)](https://travis-ci.com/kcleal/dysgu)

dysgu (pronounced *duss-gi*) is a set of tools for calling structural variants in paired-end or long read sequencing data.


Installation
------------
For convenience use the install script::

    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu
    bash INSTALL.sh

To build from source::

    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu/dysgu/htslib
    autoheader
    autoconf
    ./configure
    make
    cd ../../
    pip install -r requirements.txt
    pip install .

Run tests::

    $ dysgu test

Requires Python>=3.7, cython and >=c++11 compiler.
Python packages needed are listed in requirements.txt.

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
To call SVs, a sorted and indexed .bam/cram is needed plus an indexed reference genome in fasta format. Also a working directory must
be provided to store temporary files. There are a few ways to run dysgu depending on the type of data you have.
For paired-end data the `run` command is recommended which wraps `fetch` and `call`::

    dysgu run reference.fa temp_dir input.bam > svs.vcf

This will first call `fetch` that creates a temporary bam file and other analysis files in the working directory `samp1_temp`. These temporary files are then analysed using the `call` program.

Long reads
**********
For long-read data, the `fetch` stage may be skipped and the `call` command can be run instead - `run` is sometimes faster (PacBio Sequal II reads mainly) but involves the creation of a large
temp file::

    dysgu call --mode pacbio reference.fa temp_dir input.bam > svs.vcf
    dysgu call --mode nanopore reference.fa temp_dir input.bam > svs.vcf

The --mode=pacbio option works best with read from the SequelII platform, for older platforms use --mode=nanopore or use custom
settings.

Fetching SV reads
~~~~~~~~~~~~~~~~~
To save time, `dysgu fetch` can be run in a stream during mapping/sorting. Here, dysgu reads from stdin and
all SV associated reads will be placed in `temp_dir/temp_dir.dysgu_reads.bam`, and all input alignments will be placed in all_reads.bam::

    dysgu fetch -r all_reads.bam temp_dir -

SVs can be subsequently called using the `call` command. Additionally, the `--ibam` option is recommended for paired-end data so dysgu can infer insert
size metrics from the main alignment file. If this is not provided, dysgu will use the input.bam in the samp1_temp folder which may be less accurate::

    dysgu call --ibam all_reads.bam reference.fa temp_dir temp_dir/temp_dir.dysgu_reads.bam > svs.vcf

Alternatively, run `fetch` on an existing .bam file::

    dysgu fetch samp1_temp input.bam


Merging SVs from multiple files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiple output vcf files can be merged, e.g. tumor.vcf and normal.vcf, or illumina.vcf and pacbio.vcf::

    dysgu merge pacbio.vcf illumina.vcf > combined.vcf

For help use::

    dysgu --help
    dysgu command --help

Resource requirements
---------------------
Using a single core and depending on hard-drive speed, dysgu usually takes around 1h to analyse a 30X coverage genome of 150 bp paired-end reads and
uses < 8 GB memory. Also note that when `fetch` is utilized (run command), a large temp file is generated consisting of SV-associated reads
which can be 5 - 15 Gb in size, to remove this on completion add `--keep-temp False`.

Issues
------
If dysgu is taking a long time to run, this could be due to the complexity of the sample.
Dysgu will try and generate contigs from clusters of soft-clipped reads and remap these to the reference genome.
In this case consider increasing the `clip-length` or setting `contigs False`, or `remap False`.
Alternatively check your sample for anomalous sequences and adapter content.

If sensitivity is lower than expected for paired-end data, check that the insert size was inferred accurately, and
provide manually using the `-I` option otherwise.
