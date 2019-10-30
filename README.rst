========
dysgu-SV
========

dysgu (pronounced "*duss-key*") is a collection of tools for mapping and calling structural variants.


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

Required external tools on path: `samtools >=1.4 <http://www.htslib.org/>`_


Usage
-----
Available commands::

    $ call
    $ choose
    $ reads
    $ sv2fq
    $ view
    $ test


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

    $ dysgu call your.bam > results.csv
    $ dysgu view -f vcf results.csv > results.vcf

For help use::

    $ dysgu --help
    $ dysgu COMMAND --help

