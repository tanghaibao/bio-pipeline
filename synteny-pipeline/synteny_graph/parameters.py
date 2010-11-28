"""
Following are the parameters for the entire synteny-pipeline
"""

Tandem_Nmax = 10 # defines how "close" are tandems located (used in merge_tandems)
Nmax = 40 # maxdist between anchor points (used in synteny_score, single_linkage)
N = 3 # minimum size for cluster to report (used in single_linkage)
Cscore_cutoff = .2  # the cut-off for cscore
Synteny_cutoff = 1.5  # the cut-off for synteny_score

