========
dysgu-SV
========

dysgu (pronounced *duss-key*) is a collection of tools for calling structural variants.


Installation
------------
Install using::

    $ git clone https://github.com/kcleal/dysgu.git; cd dysgu
    $ python setup.py install

To satisfy requirements::

    $ pip install -r requirements.txt; pip install .

Run tests::

    $ dysgu test

Requires Python>=3.6, cython and >=c++11 compiler.
Python packages needed are listed in requirements.txt.


Usage
-----
Available commands::

    $ dysgu fetch   # Seperate SV reads from input bam file
    $ dysgu call    # SV calling paired-end reads
    $ dysgu view    # Merge calls from multiple dysgu call runs
    $ dysgu test    # Run basic tests


Structural variant pipeline basic usage::

    $ dysgu fetch your.bam > sv_reads.bam
    $ dysgu call sv_reads.bam > results.vcf


Fetching SV reads
~~~~~~~~~~~~~~~~~
To save time, `dysgu fetch` can be run in a stream during mapping, although here `sv_reads.bam` would
still need to be sorted prior to running `dysgu call`::

    $ bwa mem -a -t8 ref.fa read1.fq read2.fq | \
        samtools view -bh - | \
        dysgu fetch -o all_reads.bam -r sv_reads.bam  -  #  <- needs sorting

`fetch` can be run downstream of sorting which preserves sort order::

    $ samtools sort - | \
        dysgu fetch -o all_reads.bam -r sv_reads.bam  -

Alternatively, run `fetch` on an existing .bam file::

    $ dysgu fetch all_reads.bam > sv_reads.bam


Calling SVs
~~~~~~~~~~~
Input reads can be buffered which can lower run time for large memory machines. The `--buffer-size` option sets the number of alignments that will be kept in memory::

    $ dysgu call --buffer-size 10000000 sv_reads.bam > results.vcf

Calling can also be piped which can be useful for calling small regions. For this to work the `--buffer-size` must be set large enough to capture the input::

    $ samtools view -bh your.bam chr1:150000-1501000 > dysgu call - > results.vcf


Merging SVs from multiple files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If `dysgu call --out-format csv` is set, multiple output files can be merged, e.g. tumor.csv and normal.csv::

    $ dysgu view --merged-across True  tumor.csv normal.csv > pair.vcf

For help use::

    $ dysgu --help
    $ dysgu COMMAND --help

