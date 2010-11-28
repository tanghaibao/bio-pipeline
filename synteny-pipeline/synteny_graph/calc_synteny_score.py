#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
this implements a simple metric to quantify synteny score of blast hits
the score is contributed by the quality and density of the blast pairs nearby
weighted by the distance to the query blast pair

"""

import sys
import pprint
import collections
import numpy as np
from scipy.stats import norm
from scipy.spatial import cKDTree
from optparse import OptionParser
from grouper import Grouper
from parameters import Nmax


def import_genes(genes_file):

    # read in the .genes file

    ranks = {} # chromosomal position of a gene
    genes = {} # reverse of ranks, rank => gene
    print >>sys.stderr, "Import", genes_file
    fp = file(genes_file)
    for row in fp:
        chr, gene, start, stop, size, tandem_rep, label = \
                row.split()
        if gene!=tandem_rep: continue
        pos = int(label.split(":")[1])
        # bijection
        ranks[gene] = (chr, pos)
        genes[(chr, pos)] = gene
    fp.close()

    return ranks, genes


def find_nearby(hit, tree, weights, p=2):

    # the detector shape sometimes matters, in this case we use a "manhattan" detector
    dists, idxs = tree.query(hit, p=p, distance_upper_bound=Nmax, k=64)
    # I know I passed a 3-member tuple (x, y, cscore), so the distance may be off
    # but the cscore is less than one, not affecting distance (which is int'd anyway)
    neighbor_info = [(d, idx) for (d, idx) in zip(dists, idxs) if idx!=tree.n]
    # values set to tree.n are extras

    # calculate the score for each of the four quadrants
    Q1, Q2, Q3, Q4 = 0, 0, 0, 0

    for d, i in neighbor_info:
        # join the nearby hits into any member in the cluster
        neighbor = tree.data[i]
        sscore = weights[d] * neighbor[-1] # weights * cvalue
        del_x = neighbor[0]-hit[0]
        del_y = neighbor[1]-hit[1]
        if del_x==0 or del_y==0: continue # self or tandem
        if del_x > 0:
            if del_y > 0: Q1 += sscore
            else: Q4 += sscore
        else:
            if del_y > 0: Q2 += sscore
            else: Q3 += sscore

    return hit[-1] + max(Q1, Q2, Q3, Q4) # self plus neighbors


def calc_synteny_score(subfolder, species1, species2, p, weights):

    ranks1, genes1 = import_genes("data/%(species1)s.genes.filtered" % locals())
    ranks2, genes2 = import_genes("data/%(species2)s.genes.filtered" % locals())

    fp = file("%(subfolder)s/%(species1)s_%(species2)s.cscore" % locals())
    # chromosome pair => dots within a chromosome pair comparisions
    all_points = collections.defaultdict(list)
    for row in fp:
        a, b, bitscore, cscore = row.split()
        cscore = float(cscore)
        chr1, pos1 = ranks1[a]
        chr2, pos2 = ranks2[b]
        all_points[(chr1, chr2)].append((pos1, pos2, cscore))
    fp.close()

    fw = file("%(subfolder)s/%(species1)s_%(species2)s.synteny_score" % locals(), "w")
    for chr_pair in sorted(all_points.keys()):
        all_hits = all_points[chr_pair]
        
        # TODO: use different kernel shape, other than manhattan/euclidean distance
        tree = cKDTree(np.array(all_hits), leafsize=16)
        for hit in all_hits:
            synteny_score = find_nearby(hit, tree, weights, p=p)
            chr1, chr2 = chr_pair
            pos1, pos2, cscore = hit
            
            gene1 = genes1[(chr1, pos1)]
            gene2 = genes2[(chr2, pos2)]

            fw.write("%s\t%s\t" % (gene1, gene2))
            fw.write("%.3f\t%.3f\n" % (cscore, synteny_score))

    print >>sys.stderr, "synteny score written"


def normal_weights():
    # the synteny score is calculated by combining itself and adjacent blast hits
    # weighted by a function that degrades over distance
    weights = 2*(1-norm.cdf(range(Nmax+1), 0, Nmax/2))
    return weights


if __name__ == '__main__':

    weights = normal_weights()

    usage = "usage: %prog [options] species1 species2"
    parser = OptionParser(usage)
    parser.add_option("-p", "--distance_function", dest="p",
            action="store", type="int", default=2,
            help="distance function to use: 1. Manhattan dist; 2. Euclidean dist")
    
    (options, args) = parser.parse_args()

    try:
        species1 = args[0]
        species2 = args[1]
    except:
        sys.exit(parser.print_help())

    if options.p not in (1, 2):
        sys.exit(parser.print_help())

    subfolder = "data"
    calc_synteny_score(subfolder, species1, species2, options.p, weights)


