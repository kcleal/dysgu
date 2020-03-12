========
dysgu-SV
========

dysgu (pronounced *duss-key*) is a collection of tools for mapping and calling structural variants.


Installation
------------
Install using::

    $ python setup.py install
    # Or
    $ pip install -r requirements.txt; pip install .

Run tests::

    $ dysgu test

Requires Python>=3.6, cython and >=c++11 compiler.
Python packages needed are listed in requirements.txt.


Usage
-----
Available commands::

    $ dysgu call    # SV calling paired-end reads
    $ dysgu choose  # Choose paired-end or single read alignments from sam input
    $ dysgu fetch   # Seperate SV reads from input bam file
    $ dysgu view    # Merge calls from multiple dysgu call runs
    $ dysgu test    # Run basic tests


Choosing alignments from a list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Given a list of candidate alignments for paired-end reads, `dysgu choose` will perform split read-pairing, or for single
end reads will selects an optimal spanning set of alignments.

For pairing using paired-end alignments from a candidate list::

    $ dysgu choose all_alignments.sam > output_alignments.sam

Or can be run in a stream using bwa mem in all-mapping mode (use of -a option in bwa)::

    $ bwa mem -a -t8 ref.fa read1.fq read2.fq | dysgu choose - > output_alignments.sam

Or run in single end mode::

    $ bwa mem -a -t8 ref.fa contigs.fa | dysgu choose --paired False - > output_alignments.sam


Calling structural variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Structural variant pipeline basic usage::

    $ dysgu fetch your.bam > sv_reads.bam
    $ dysgu call sv_reads.bam > results.vcf

To save time, `dysgu fetch` can be piped during mapping::

    $ bwa mem -a -t8 ref.fa read1.fq read2.fq | \
        samtools view -bh - | \
        dysgu fetch -o all_reads.bam -r sv_reads.bam
    $ dysgu call sv_reads.bam > results.vcf

Input reads can also be buffered which can lower run time for large memory machines. The `--buffer-size` option sets the number of alignments that will be kept in memory::

    $ dysgu call --buffer-size 10000000 sv_reads.bam > results.vcf

Calling can also be piped which can be useful for calling small regions. For this to work the `--buffer-size` must be set large enough to capture the input::

    $ samtools view -bh your.bam chr1:150000-1501000 > dysgu call - > results.vcf

If `dysgu call --out-format csv` is set, multiple output files can be merged, e.g. tumor.csv and normal.csv::

    $ dysgu view --merged-across True  tumor.csv normal.csv > pair.vcf

For help use::

    $ dysgu --help
    $ dysgu COMMAND --help

