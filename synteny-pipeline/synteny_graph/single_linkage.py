#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Use single linkage to find out synteny clusters
# We need to set two parameters
# Nmax - max distance between any two hits to be chained
# N - minimum size of the cluster to report

import sys
import pprint
import collections
import operator
import numpy as np
from optparse import OptionParser
from grouper import Grouper
from parameters import Nmax, N, Cscore_cutoff

def import_genes(filename):

    # read in the .genes file

    ranks = {} # chromosomal position of a gene
    print >>sys.stderr, "Import", filename
    gene_set = set()
    fp = file(filename)
    for row in fp:
        chr, gene, start, stop, size, tandem_rep, label = \
                row.split()
        label = int(label.split(":")[1])
        ranks[gene] = (chr, label)
    return ranks

def single_linkage(points, max_dist=Nmax, min_cluster_size=N):
    """
    points are (x-index, y-index, cscore) per chromosome pair.
    """
    # This is the core single linkage algorithm
    # this behaves in O(n) complexity: we iterate through the pairs, for each pair
    # we look back on the adjacent pairs to find links

    clusters = Grouper()
    n = len(points)
    points.sort()
    for i in xrange(n):
        for j in xrange(i-1, -1, -1):
            # x-axis distance
            del_x = points[i][0]-points[j][0]
            if del_x > max_dist: break
            # y-axis distance
            del_y = points[i][1]-points[j][1]
            if del_x + abs(del_y) > max_dist: continue
            #if abs(del_y) > Nmax: continue
            # otherwise join
            clusters.join(points[i], points[j])
    clusters = [cluster for cluster in list(clusters) if len(cluster)>=min_cluster_size]
    return clusters


def write_clusters(filehandle, clusters, chr_pair):
    for cluster in clusters:
        cluster_score = sum(x[2] for x in cluster)
        filehandle.write("# cluster score %.3f \n" % (cluster_score)) 
        # these are seed anchors (mutual best)
        for pos_pair in cluster:
            chr1, chr2 = chr_pair
            pos1, pos2, synteny_score = pos_pair
            filehandle.write("%s\t%d\t" % (chr1, pos1) )
            filehandle.write("%s\t%d\t%.3f\n" % (chr2, pos2, synteny_score) )


def read_clusters(filename):
    fp = file(filename)
    row = fp.readline()
    cluster_id = 0
    clusters = []
    while row:
        lines = []
        row = fp.readline()
        while row and row[0]!="#":
            chr1, pos1, chr2, pos2, cscore = row.split()
            pos1, pos2 = int(pos1), int(pos2)
            cscore = float(cscore)
            lines.append(((chr1, pos1), (chr2, pos2), cscore))
            row = fp.readline()
        lines.sort() # this will sort with respect to the first axis
        clusters.append(lines)
    return clusters


def distance(gene1, gene2):
    chr1, pos1 = gene1
    chr2, pos2 = gene2
    if chr1 != chr2: return Nmax + 1 # this ensures un-chainable
    return abs(pos1 - pos2)


def distance_x(cluster_i, cluster_j):
    # x-sorted
    min_i_x, max_i_x = min(cluster_i)[0], max(cluster_i)[0]
    min_j_x, max_j_x = min(cluster_j)[0], max(cluster_j)[0]

    del_x1 = distance(min_i_x, max_j_x)
    del_x2 = distance(max_i_x, min_j_x)
    del_x3 = distance(min_i_x, min_j_x)
    del_x4 = distance(max_i_x, max_j_x)
    
    return min(del_x1, del_x2, del_x3, del_x4)

def distance_y(cluster_i, cluster_j):
    # y-sorted
    min_i_y = min(cluster_i, key=operator.itemgetter(1))[1]
    max_i_y = max(cluster_i, key=operator.itemgetter(1))[1]
    min_j_y = min(cluster_j, key=operator.itemgetter(1))[1]
    max_j_y = max(cluster_j, key=operator.itemgetter(1))[1]

    del_y1 = distance(min_i_y, max_j_y)
    del_y2 = distance(max_i_y, min_j_y)
    del_y3 = distance(min_i_y, min_j_y)
    del_y4 = distance(max_i_y, max_j_y)

    return min(del_y1, del_y2, del_y3, del_y4)

def merge_clusters(chain, clusters):

    # there are, in general, two kinds of breakpoints
    # those that are induced by inversions, and those by translocations
    # inversion-breakpoints are excessive breakpoints that I want to remove
    
    chain_num = len(chain)
    mergeables = Grouper() # disjoint sets of clusters that can be merged
    for j in xrange(chain_num):
        cj = chain[j]
        mergeables.join(cj, cj)
        for i in xrange(j-1, -1, -1):
            ci = chain[i]
            del_x = distance_x(clusters[ci], clusters[cj])
            if del_x > Nmax: continue 

            del_y = distance_y(clusters[ci], clusters[cj])
            if del_x + del_y > Nmax: continue
            mergeables.join(ci, cj)

    to_merge = {} 
    for mergeable in mergeables:
        for m in mergeable:
            to_merge[m] = min(mergeables[m])

    merged_chain = []
    for c in chain:
        if to_merge[c]==c: # i.e. parent of mergeables
            merged_chain.append(c)

    # refresh clusters list, merge chains
    for k, v in to_merge.iteritems():
        if to_merge[k]!=k: # i.e. not map to self
            clusters[v].extend(clusters[k])

    # maintain the x-sort
    [cluster.sort() for cluster in clusters]

    # nothing is merged
    updated = (len(merged_chain) != chain_num)
    return merged_chain, updated


def recursive_merge_clusters(chain, clusters):

    # as some rearrangment patterns are recursive, the extension of blocks
    # will take several iterations
    while 1: 
        chain, updated = merge_clusters(chain, clusters)
        print >>sys.stderr, "merging..."
        if not updated: break

    return chain, clusters


def write_chain(fw_chain, chain, clusters):

    # output the anchor points in each chain
    for c in chain:
        lines = clusters[c]
        print >>fw_chain, "# mergedcluster score %.3f" % sum(x[-1] for x in lines)
        for line in lines:
            gene1, gene2, synteny_score = line
            chr1, pos1 = gene1
            chr2, pos2 = gene2
            print >>fw_chain, \
                    "%s\t%d\t%s\t%d\t%.3f" % (chr1, pos1, chr2, pos2, synteny_score)


if __name__ == '__main__':
        
    usage = "usage: %prog [options] species1 species2 pairs_file"
    parser = OptionParser(usage)
    parser.add_option("-m", "--mergeblocks", action="store_true",
            dest="mergeblocks", default=False, 
            help="merge adjacent blocks after considering inversions")

    (options, args) = parser.parse_args()

    try:
        species1 = args[0]
        species2 = args[1]
        pairs_file = args[2]
    except:
        sys.exit(parser.print_help())

    self_match = (species1==species2)
    print >>sys.stderr, "comparing self"

    ranks1 = import_genes("data/%(species1)s.genes.filtered" % locals())
    if self_match:
        ranks2 = ranks1
    else:
        ranks2 = import_genes("data/%(species2)s.genes.filtered" % locals())

    fp = file(pairs_file)
    # chromosome pair => dots within a chromosome pair comparisions
    all_points = collections.defaultdict(list)
    good_points = collections.defaultdict(list)
    for row in fp:
        
        atoms = row.split()
        a, b = atoms[:2]
        synteny_score = atoms[-1]
        if synteny_score=="n.a.": 
            synteny_score = Cscore_cutoff
        else:
            synteny_score = float(synteny_score)
        if synteny_score < 1.5: continue

        chr_a, label_a = ranks1[a]
        chr_b, label_b = ranks2[b]

        # if comparing self, keep only the ones in lower triangle
        if self_match and (chr_a, label_a) > (chr_b, label_b): continue

        good_points[(chr_a, chr_b)].append((label_a, label_b, synteny_score))

    fp.close()

    fw = file("data/%(species1)s_%(species2)s.cluster" % locals(), "w")
    for chr_pair, points in good_points.iteritems():
        clusters = single_linkage(points)
        write_clusters(fw, clusters, chr_pair)

    fw.close()
    print >>sys.stderr, "clusters found and written to .cluster file"

    # additional efforts to merge rearranged blocks (reduce breakpoints)
    if not options.mergeblocks: sys.exit(0) 

    # `chain' only contains integer ID that refer to the index in the clusters
    # `clusters' contain the actual gene pairs

    clusters = read_clusters("data/%(species1)s_%(species2)s.cluster" % locals)
    chain = range(len(clusters)) # the original chain contains everything
    chain, clusters = recursive_merge_clusters(chain, clusters)
    fw_chain = file("data/%(species1)s_%(species2)s.mergedcluster", "w")
    write_chain(fw_chain, chain, clusters)

    print >>sys.stderr, "mergedclusters found and written to .mergedcluster file"

