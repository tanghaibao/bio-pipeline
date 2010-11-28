Description
--------------
This package is only planned while not yet usable, with a list of TO-DOs are scheduled, and help document will be updated with development

- filter the pairwise mafs to get non-overlapping set, based on a reference sequence
- implement quota-alignment algorithm to deal with the mafs
- combine pairwise mafs to construct multiple mafs


Scripts
-----------
Caution: none of the following scripts are generic enough to be re-used.

- ``maf_to_bed.py``, converts the ``.maf`` file to the ``.bed`` format for downstream analysis (like coverage, or intersection with annotation features)
- ``maf_chains.py``, filters the pairwise mafs to get two 2-monotonous chains
- ``maf_select.py``, select a pair of chromosomes and sends all alignments to file
- ``maf_to_blast8.py``, convert the .maf formatted file to the BLAST m8 format
- ``axt_chain/run.sh``, bash script to use the UCSC BLASTZ/CHAIN/NET pipeline, by following instructions at `here <http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto>`_.

