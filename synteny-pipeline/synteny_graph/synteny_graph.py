#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""

Consider a graph, vertices are individual genes from different species, 
and edges are synteny score, reconstruct multi-gene families that match certain 
expected ratios (usually expected copy numbers following polyploidy)

So the input is a bunch of pairwise blast hits already filtered by synteny, the
algorithm follows:

- sort all the pairs in decreasing strength of synteny
- merging genes into gene family with some constraints:
    1. respect quota set up by the polyploidies
    2. very close genes are only counted once in the family \
            (might represent split gene models or long distance local dup)

"""

import os
import sys
import pprint
import collections
from glob import glob
from optparse import OptionParser
from grouper import Grouper
from parameters import Tandem_Nmax, Synteny_cutoff

from scipy.stats import norm
from scipy.spatial import cKDTree
from calc_synteny_score import calc_synteny_score, \
        normal_weights


def read_genes(genes_file, species, all_ranks):

    print >>sys.stderr, "Import", genes_file 
    fp = file(genes_file)
    for row in fp:
        chr, gene, start, stop, size, tandem_rep, label = \
                row.split()
        pos = int(label.split(":")[1])
        all_ranks[gene] = (species, chr, pos)
    fp.close()


def read_cscore(score_file):

    print >>sys.stderr, "Import", score_file 

    all_cscore = {}
    fp = file(score_file)
    for row in fp:
        gene_a, gene_b, bitscore, cscore = row.split() 
        cscore = float(cscore)
        all_cscore[(gene_a, gene_b)] = row

    return all_cscore


def filter_cscore(family_file, score_file, filtered_file):

    all_cscore = read_cscore(score_file)

    fp = file(family_file)
    fw = file(filtered_file, "w")

    print >>sys.stderr, "Make refined scores %s" % filtered_file
    for row in fp:
        genes = row.strip().split(",")
        gene_nums = len(genes)
        for i in xrange(gene_nums):
            genei = genes[i]
            for j in xrange(i+1, gene_nums):
                genej = genes[j]
                if (genei, genej) in all_cscore:
                    fw.write(all_cscore[(genei, genej)])
                if (genej, genei) in all_cscore:
                    fw.write(all_cscore[(genej, genei)])

    fp.close()
    fw.close()


def make_new_scores(subfolder, species_list, family_file):

    weights = normal_weights()
    p = 2

    ngenomes = len(species_list) # number of genomes
    for i in xrange(ngenomes):

        # if the quota asks for duplication
        speciesi = species_list[i]
        if quota[speciesi] > 1:
            score_file = "data/%s_%s.cscore"%(speciesi, speciesi)
            filtered_file = score_file.replace("data", subfolder)
            filter_cscore(family_file, score_file, filtered_file)
            calc_synteny_score(subfolder, speciesi, speciesi, p, weights)

        for j in xrange(i+1, ngenomes):
            speciesj = species_list[j]
            score_file = "data/%s_%s.cscore"%(speciesi, speciesj)
            filtered_file = score_file.replace("data", subfolder)
            filter_cscore(family_file, score_file, filtered_file)
            calc_synteny_score(subfolder, speciesi, speciesj, p, weights)
    

def read_synteny_score(score_file, gene_pairs):

    print >>sys.stderr, "Import", score_file 

    fp = file(score_file)
    for row in fp:
        gene_a, gene_b, cscore, synteny_score = row.split() 
        synteny_score = float(synteny_score)
        if synteny_score < Synteny_cutoff: continue
        gene_pairs.append((synteny_score, gene_a, gene_b))


def mergeable(group1, group2, all_ranks, quota):

    # rule no.1- respect quota
    # rule no.2- but only count close genes once
    micro_grouper = Grouper() # data structure to check rule no.2
    merged_group = group1 + group2 # attempted merge
    nmerged = len(merged_group)

    # do all pairwise comparisons to find closely located genes
    # TODO: silly implementation, not efficient
    for i, genei in enumerate(merged_group):
        speciesi, chri, posi = all_ranks[genei]
        micro_grouper.join(genei)
        for j in xrange(i+1, nmerged):
            genej = merged_group[j]
            speciesj, chrj, posj = all_ranks[genej]
            if speciesi==speciesj and chri==chrj and abs(posi-posj)<=Tandem_Nmax/2:
                micro_grouper.join(genei, genej)

    species_count = collections.defaultdict(int) # data structure to check rule no.1
    for gene_group in micro_grouper:
        species = all_ranks[gene_group[0]][0]
        species_count[species] += 1

    for species, count in species_count.items():
        if count>quota[species]: 
            return False

    return True


def make_family(gene_pairs, all_ranks, quota):
    
    print >>sys.stderr, "... gene family clustering started"

    g = Grouper() 

    gene_pairs.sort(reverse=True)
    #pprint.pprint(gene_pairs[:10])
    for synteny_score, gene1, gene2 in gene_pairs:
        # attempt to join the two genes

        g.join(gene1)
        g.join(gene2)
        group1, group2 = g[gene1], g[gene2]

        if mergeable(group1, group2, all_ranks, quota):
            g.join(gene1, gene2)

    return g


def write_family(family_file, gene_family):

    fw = file(family_file, "w")
    sorted_family = [] # for better output
    for family in gene_family:
        if len(family)==1: continue # do not write singletons
        sorted_family.append(sorted(family))
    sorted_family.sort()

    for family in sorted_family:
        fw.write("%s\n" % (",".join(family)))
    print >>sys.stderr, len(sorted_family), "gene families written to ", family_file 


def synteny_graph(subfolder, species_list, quota):
    all_ranks = {} # genename => position
    # read in all the position information
    for species in species_list:
        read_genes("data/%s.genes.filtered" % species, species, all_ranks)

    ngenomes = len(species_list) # number of genomes
    # read in all the syntenic gene pairs
    gene_pairs = []
    for i in xrange(ngenomes):

        # if the quota asks for duplication
        speciesi = species_list[i]
        if quota[speciesi] > 1:
            score_file = "%s/%s_%s.synteny_score"%(subfolder, speciesi, speciesi)
            read_synteny_score(score_file, gene_pairs)

        for j in xrange(i+1, ngenomes):
            speciesj = species_list[j]
            score_file = "%s/%s_%s.synteny_score"%(subfolder, speciesi, speciesj)
            read_synteny_score(score_file, gene_pairs)

    print >>sys.stderr, len(gene_pairs), "syntenic pairs imported"
    
    g = make_family(gene_pairs, all_ranks, quota)
    return g


if __name__ == '__main__':

    usage = "usage: %prog [options] -s species_list"
    parser = OptionParser(usage)
    parser.add_option("-s", "--species", dest="species",
            action="store", type="string", default="all", 
            help="build gene family from a list of species, separated by ':' "\
                    "(like athaliana:grape), default is to all in data/ folder "\
                    "this needs to be sorted, so athaliana goes before grape")
    parser.add_option("-q", "--quota", dest="quota",
            action="store", type="string", default="4:1:1:2", 
            help="the expected gene copies following the common ancestor "\
                    "(for example athaliana:grape, should be 4:1),"\
                    "this has to match the order you specified in the -s !")
    parser.add_option("-r", "--refine", dest="refine", 
            action="store_true", default=False,
            help="further refine the synteny graph to reduce false positives "\
                    "(will calculate synteny score recursively)")

    (options, args) = parser.parse_args()

    if options.species=="all":
        species_list = [os.path.basename(x).split(".")[0] for x in \
                glob("data/*.genes.filtered")]
    else:
        try:
            species_list = options.species.split(":")
        except:
            sys.exit(parser.print_help())

    species_list.sort()
    try:
        quota_list = [int(x) for x in options.quota.split(":")]
    except:
        sys.exit(parser.print_help())

    assert len(species_list)==len(quota_list), "number of genomes given do not match!"

    quota = dict(zip(species_list, quota_list))

    for s in sorted(quota.keys()): 
        print >>sys.stderr, s, quota[s], "copies"

    subfolder = "data"
    g = synteny_graph(subfolder, species_list, quota)
    family_file = "rosid_quota_family"
    write_family(family_file, g)

    if not options.refine: sys.exit(0)

    # Further refinement --
    # Here is the idea, dis-assemble the family yet again
    # run cscore->synteny_score->construct_graph multiple times
    # noise will get away in a few iterations

    family_size_unrefined = len(g)+1
    for j in xrange(1, 100):

        print >>sys.stderr, "=== Refinement iteration %d ===" % j
        family_size_unrefined = len(g)
        subfolder = "work"
        make_new_scores(subfolder, species_list, family_file)
        g = synteny_graph(subfolder, species_list, quota)
        if len(g) >= family_size_unrefined: # convergence
            break
        write_family(family_file, g)


