Synteny graph is consisted of several utility scripts that can be run sequentially through a batch shell script.

Utility scripts
===============
- ``merge_tandem.py`` looks for the self blast and then produce a gene order with the tandems merged
- ``calc_cscore.py`` calculates the cscore for the blast results
- ``calc_synteny_score.py`` calculates the synteny score for the blast results
- ``synteny_graph.py`` clusters the genes within the filtered syntenic blocks and consolidates into a gene group

