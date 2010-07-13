Ks calculation script
======================

Installation
-------------
The python script will run on Linux operation system. To install, you must have the following softwares in place, please contact system admin if you do not know how to install,

* `biopython <http://www.biopython.org>`_
* `CLUSTALW <http://www.ebi.ac.uk/Tools/clustalw2/index.html>`_
* `PAL2NAL <http://coot.embl.de/pal2nal/>`_
* `PAML <http://abacus.gene.ucl.ac.uk/software/paml.html>`_

Please remember the installation path for CLUSTALW, PAL2NAL and PAML. You will then need to modify the script ``synonymous_calc.py`` line ``19-21`` to change the path.

Preparing data
---------------
You need to prepare two FASTA files, one file with protein seqs, one file with gene seqs, they contain the gene pairs::

    >gene1
    ATAGATATATATA
    >gene2
    ATATAGAGAGAGA
    >gene3
    AGAGAGAGAGAGA
    >gene4
    ATAGAGAGAGAGA

This will calculate two pairs: gene1-gene2 and gene3-gene4. Make sure that your protein seq file corresponds to your gene seq file, in exactly the same order.

Usage
------
Finally, run the command like this::

    $ python synonymous_calc.py test.pep test.cds test.ks

where ``test.pep`` is your protein file, ``test.cds`` is your CDS file and your result is in ``test.ks``. The result is a comma-delimited file, you can open it in EXCEL, columns correspond to::

    Pair_ID; Yang-Nielson method Ks, Yang-Nielson method Ka, Nei-Gojobori method Ks, Nei-Gojobori method Ka

I personally recommend Nei-Gojobori method, from past experience. `Nei-Gojobori method <http://www.megasoftware.net/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_modified_nei_gojobori_method.htm>`_ is a simple correction for multiple substitutions.
