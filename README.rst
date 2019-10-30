====
fnfi
====

fnfi is a collection of tools for mapping and calling structural variants.


Installation
------------
Install using::

    $ python setup.py install
    # Or
    $ pip install .

Requires Python >= 3, c++ compatible compiler. Python packages needed: click,
numpy, pandas, pysam, pybedtools, natsort, networkx, scikit-learn, ncls.

Required external tools on path: `samtools >=1.4 <http://www.htslib.org/>`_

Recommended external tools on path: `bwa mem <https://github.com/lh3/bwa/>`_

Usage
-----
Available commands::

    $ fnfi run
    $ fnfi find-reads
    $ fnfi align
    $ fnfi call-events


Choosing alignments from a list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Given a list of candidate alignments for paired-end reads, `fnfi align` will perform split read-pairing, or for single
end reads will selects an optimal spanning set of alignments.

For pairing using paired-end alignments from a candidate list::


    $ fnfi align all_alignments.sam > output_alignments.sam

Or can be run in a stream using bwa mem in all-mapping mode (use of -a option in bwa)::


    $ bwa mem -a -t8 ref.fa read1.fq read2.fq | fnfi align - > output_alignments.sam

Or run in single end mode::


    $ bwa mem -a -t8 ref.fa contigs.fa | fnfi align --paired False - > output_alignments.sam


Calling structural variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Currently under development and not intended for use.
Structural variant pipeline basic usage::

    $ fnfi run reference.fa your.bam > results.csv

For help use::

    $ fnfi --help
    $ fnfi COMMAND --help

