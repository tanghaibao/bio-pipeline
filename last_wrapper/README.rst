How to use
===========

$ python last.py --path=. -a 32 lotus.cds lotus.cds -o lotus.lotus.last

The --path should contain the PATCHED `lastal` that is capable to generate BLAST
tabular output. The patch is also shown in the folder, fix the lastal source
code accordingly if you need to recompile.

Drop scripts `last.py` and `last-helper.py` (functions used by last.py) anywhere
you want but perhaps you want to put it in quota-align folder. Drop
`lastal/lastdb/lastex` anywhere you want but specify the path with --path.

Biopython may have been installed already, but I copied it over just in case.
