Python GFF sequence loader
===========================

Setup
------
This python script depends on external libraries, namely:

* `GFFutils <http://github.com/daler/GFFutils>`__, for parsing GFF file
* `pyfasta <http://github.com/brentp/pyfasta>`__, for extracting subsequence

Both packages are easy-installable::

    $ easy_install GFFutils
    $ easy_install pyfasta

Usage
------
Just locate the GFF file and the FASTA file, and then run the script::

    $ ./gff_loader.py athaliana.gff athaliana.fa 

This by default, extracts the ``mRNA`` IDs and pull out all subfeatures (``CDS``)
that have the same parent and concatenate the seqs together, reverse-complement
if needed. However sometimes the GFF does not use the standard names, this is
when you need to do::

    $ ./gff_loader.py athaliana.gff athaliana.fa --parents Gene --children exon

If the GFF has meant coding sequences, but uses a different term.
