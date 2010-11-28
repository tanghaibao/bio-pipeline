#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Take self-blast results, remove all tandems that can be connected 
within a certain specified distances

Rank all the genes after tandem removal, and note which tandem 
group it belongs (represented by a longest peptide in the group)

"""

import sys
import collections
from grouper import Grouper, gene_name
from parameters import Tandem_Nmax
from optparse import OptionParser

def load_sizes(fp_sizes):

    # load the gene size info (for keeping longest gene in a tandem group)

    fp_sizes.seek(0)
    sizes = {} # gene => size
    print >>sys.stderr, "Read .sizes file"
    for row in fp_sizes:
        gene, size = row.split()
        # sizes are calculated for the transcripts, and we keep longest transcript
        gene, size = gene_name(gene), int(size)
        if gene not in sizes or size > sizes[gene]:
            sizes[gene] = size
    return sizes

def load_geneorders(fp_gff):

    # load gene orders before any filtering
    
    fp_gff.seek(0)
    tandem = Grouper()
    print >>sys.stderr, "Read .genes file"
    # chromosome => gene_list in that chromosome
    chr_ranks = collections.defaultdict(list)
    ranks = {} # gene => rank postion
    for row in fp_gff:
        chr, gene, start, stop = row.split()
        start = int(start)
        chr_ranks[chr].append((start, gene, chr))
        tandem.join(gene)
    for v in chr_ranks.itervalues():
        gene_rank = 0
        for start, gene, chr in sorted(v): 
            ranks[gene] = (chr, gene_rank)
            gene_rank += 1
    return ranks, tandem

def process_tandems(fp_blast, sizes, ranks, tandem):

    fp_blast.seek(0)
    total_lines = sum(1 for row in fp_blast)
    print >>sys.stderr, \
            "Read self BLAST file (total %d lines)" % total_lines
    fp_blast.seek(0)
    j = 0
    for row in fp_blast:
        
        j+=1
        if j % 100000 == 0: print >>sys.stderr, j, "read..."

        atoms = row.split()
        a, b = atoms[:2]
        a, b = gene_name(a), gene_name(b)
        if a not in ranks or b not in ranks: continue
        chr_a, rank_a = ranks[a]
        chr_b, rank_b = ranks[b]
        if chr_a==chr_b and abs(rank_a-rank_b) <= Tandem_Nmax: 
            tandem.join(a, b)

    tandem_removed = set() # the filtered gene set
    tandem_map = {} # unfiltered => filtered
    for tandem_group in tandem:
        longest_gene, longest_size = "", 0
        for gene in tandem_group:
            if gene in sizes:
                gene_size = sizes[gene]
                if gene_size > longest_size:
                    longest_gene, longest_size = gene, gene_size 
        for gene in tandem_group:
            tandem_map[gene] = longest_gene
        tandem_removed.add(longest_gene)

    print >>sys.stderr, len(tandem_removed), "genes after tandem removal"

    return tandem_map, tandem_removed

def reorder(fp_gff, sizes, tandem_removed):

    # re-order everything after removing tandems and non-protein genes

    fp_gff.seek(0)
    _chr = ""
    pos_map = {}
    for row in fp_gff:
        chr, gene, start, stop = row.split()
        # mostly if gene is not in the sizes (fasta), it is ncRNA
        if gene not in sizes: continue

        if chr!=_chr: 
            j = 0
        if gene in tandem_removed:
            new_pos = "%s:%05d"%(chr, j)
            pos_map[gene] = new_pos
            j += 1 # increment only for the non-tandem genes
        _chr = chr

    return pos_map

def write_neworder(fp_gff, fw_filtered, sizes, tandem_map, pos_map):
    # write the new order and tandem representative to a file
    print >>sys.stderr, "Write filtered .genes files"
    fp_gff.seek(0)
    for row in fp_gff:
        chr, gene, start, stop = row.split()
        # mostly if gene is not in the fasta, it is miRNA
        if gene not in sizes: continue

        tandem_rep = tandem_map[gene]
        fw_filtered.write("\t".join([row.strip(), \
                str(sizes[gene]), tandem_rep, pos_map[tandem_rep]])) 
        fw_filtered.write("\n")


if __name__ == '__main__':
        
    usage = "usage: %prog [options] species"
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()

    try:
        species = args[0] # species name
    except:
        parser.print_help()
        sys.exit(1)

    annotation_folder = "/home/bao/data/annotations"
    blast_folder = "/home/bao/blast/results"

    fp_sizes = file("%(annotation_folder)s/%(species)s.pep.sizes" % locals())
    fp_gff = file("%(annotation_folder)s/%(species)s.genes" % locals())
    fp_blast = file("%(blast_folder)s/%(species)s_%(species)s.blastz" % locals()) 
    fw_filtered = file("data/%(species)s.cds.filtered" % locals(), "w")

    sizes = load_sizes(fp_sizes)
    ranks, tandem = load_geneorders(fp_gff)
    
    tandem_map, tandem_removed = process_tandems(fp_blast, sizes, ranks, tandem)
    pos_map = reorder(fp_gff, sizes, tandem_removed)
    
    write_neworder(fp_gff, fw_filtered, sizes, tandem_map, pos_map)

    # close all files after finished
    fp_sizes.close()
    fp_gff.close()
    fp_blast.close()
    fw_filtered.close()
