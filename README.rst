========
dysgu-SV
========

dysgu (pronounced *duss-key*) is a collection of tools for calling structural variants.


Installation
------------
For convenience use the install script::

    git clone --recursive https://github.com/kcleal/dysgu.git
    cd dysgu
    bash INSTALL.sh -j4

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

Requires Python>=3.6, cython and >=c++11 compiler.
Python packages needed are listed in requirements.txt.


Usage
-----
Available commands::

    dysgu run     # Run using default arguments, wraps fetch and call commands
    dysgu fetch   # Seperate SV reads from input bam file
    dysgu call    # SV calling
    dysgu view    # Merge calls from multiple dysgu call runs
    dysgu test    # Run basic tests

Structural variant pipeline basic usage::

    dysgu run your.bam > calls.vcf


Fetching SV reads
~~~~~~~~~~~~~~~~~
To save time, `dysgu fetch` can be run in a stream during mapping, although here `sv_reads.bam` would
still need to be sorted prior to running `dysgu call`::

    bwa mem -a -t8 ref.fa read1.fq read2.fq | \
        samtools view -bh - | \
        dysgu fetch -r all_reads.bam -o sv_reads.bam  -  #  -r and -o still need sorting

`fetch` can be run downstream of sorting which preserves sort order::

    samtools sort your.bam | \
        dysgu fetch -r your.bam -o sv_reads.bam  -

Alternatively, run `fetch` on an existing .bam file::

    dysgu fetch your.bam > sv_reads.bam


Calling SVs
~~~~~~~~~~~
For a sorted .bam or .cram file, use the `run` command with -p/--procs to specify reading threads::

    dysgu run -p 8 your.bam > svs.vcf

Input reads can be buffered which can lower run time for large memory machines. The `--buffer-size` option sets the number of alignments that will be kept in memory::

    dysgu call --buffer-size 10000000 sv_reads.bam > results.vcf

Calling can also be piped which can be useful for calling small regions. For this to work the `--buffer-size` must be set large enough to capture the input::

    samtools view -bh your.bam chr1:150000-1501000 > dysgu call - > results.vcf


Merging SVs from multiple files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiple output files can be merged, e.g. tumor.vcf and normal.vcf, or illumina.vcf and pacbio.vcf::

    dysgu view tumor.csv normal.csv > pair.vcf

For help use::

    dysgu --help
    dysgu command --help

