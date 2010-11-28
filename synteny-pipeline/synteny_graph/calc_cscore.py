#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use c-value filtering
see sea anemone and amphioxus genome paper;

c-value is a float number in (0, 1], defined for each blast pair

c-value = score(a, b) / max(score(a, x), score(y, b))
for any x in genome X and any y in genome Y;

in other words, reciprocal best would have c-value of 1
anything lower than 1 generalize the concept of reciprocal best,
but still requires a pretty high score 
"""

import sys
import operator
from optparse import OptionParser
from grouper import Grouper, gene_name
from parameters import Cscore_cutoff 

class Hit:

    # Keep track of the value for the hit in the blast pairs
    # hit gene => bitscore

    def __init__(self, name, score):
        self.mapping = {}
        self.mapping[name] = score

    def __len__(self):
        return len(self.mapping)

    def update(self, name, score):
        # Keep the best score for the current hit
        if name not in self.mapping or \
                score > self.mapping[name]: 
            self.mapping[name] = score

    def calc_score(self):
        self.best_score = max(self.mapping.values())


def import_genes(filename, tandem_map):

    print >>sys.stderr, "Import", filename
    gene_set = set()
    fp = file(filename)
    for row in fp:
        chr, gene, start, stop, size, tandem_rep, label = \
                row.split()
        tandem_map[gene] = tandem_rep
        if gene!=tandem_rep: continue
        gene_set.add(tandem_rep)
    return gene_set


def import_blast(filename):

    fp = file(filename)
    total_lines = sum(1 for row in fp)
    print >>sys.stderr, \
            "Read BLAST file %(filename)s (total %(total_lines)d lines)" % locals()
    fp.seek(0)
    j = 0
    for row in fp:

        j += 1
        if j % 100000 == 0: print >>sys.stderr, j, "read..."
            
        atoms = row.split()
        a, b, bitscore = atoms[0], atoms[1], float(atoms[-1])
        a, b = gene_name(a), gene_name(b)
        if a not in tandem_map or b not in tandem_map:
            continue
        a, b = tandem_map[a], tandem_map[b]

        if a == b: continue
        # keep the best blast hit
        if a not in blast_pool:
            blast_pool[a] = Hit(b, bitscore)
        else:
            blast_pool[a].update(b, bitscore)

    fp.close()


if __name__ == '__main__':
        
    usage = "usage: %prog [options] species1 species2"
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()

    try:
        species1 = args[0]
        species2 = args[1]
    except:
        sys.exit(parser.print_help())

    blast_folder = "/home/bao/blast/results"
    genes_folder = "data"

    blast_pool = {} 
    tandem_map = {}

    print >>sys.stderr, "Process all blast hits for clustering"
    gene_set1 = import_genes("%(genes_folder)s/%(species1)s.genes.filtered" % \
            locals(), tandem_map)
    gene_set2 = import_genes("%(genes_folder)s/%(species2)s.genes.filtered" % \
            locals(), tandem_map)

    cross_blast1 = "%(blast_folder)s/%(species1)s_%(species2)s.blastp" % locals()
    cross_blast2 = "%(blast_folder)s/%(species2)s_%(species1)s.blastp" % locals()

    import_blast(cross_blast1)
    import_blast(cross_blast2)

    print >>sys.stderr, len(blast_pool), "records read"
    [v.calc_score() for v in blast_pool.itervalues()]

    print >>sys.stderr, "Filter BLAST file"
    fp = file(cross_blast1)
    fw = file("data/%(species1)s_%(species2)s.cscore" % locals(), "w")
    for a, hit in sorted(blast_pool.items()):
        for b, score in hit.mapping.items():
            if a not in blast_pool or b not in blast_pool: continue
            if a not in gene_set1: continue
            
            # see __doc__ for explanation of this formula
            #c_score = score/blast_pool[a].best_score  # one direction only 
            c_score = score/max([blast_pool[a].best_score, blast_pool[b].best_score])
            
            # c_score filter here
            if c_score < Cscore_cutoff: continue

            fw.write("%(a)s\t%(b)s\t%(score).1f\t%(c_score).3f\n" % locals() )

    fp.close()
    fw.close()
