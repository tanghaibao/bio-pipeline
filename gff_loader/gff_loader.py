#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
%prog gff_file fasta_file [--options]

Parses the selected features out of GFF, with subfeatures concatenated together.
For example, to get the CDS sequences, do this::
    %prog athaliana.gff athaliana.fa --parents mRNA --children CDS
'''

import sys
import os
import os.path as op
import itertools
import GFFutils

from pyfasta import Fasta

def main(gff_file, fasta_file, parents, children):

    db_file = gff_file + ".db"

    if not op.exists(db_file):
        GFFutils.create_gffdb(gff_file, db_file)

    f = Fasta(fasta_file)
    g = GFFutils.GFFDB(db_file)

    parents = set(parents.split(','))
    parents_iter = [g.features_of_type(x) for x in parents]
    parents_list = itertools.chain(*parents_iter)
    children_list = set(children.split(','))

    for feat in parents_list:

        children = []
        for c in g.children(feat.id, 1):

            if c.featuretype not in children_list: continue
            child = f.sequence(dict(chr=c.chrom, start=c.start, stop=c.stop,
                strand=c.strand))
            children.append((child, c))

        if not children: 
            print >>sys.stderr, "[warning] %s has no children with type %s" \
                                    % (feat.id, ','.join(children_list))
            continue
        # sort children in incremental position
        children.sort(key=lambda x: x[1].start)
        # reverse children if negative strand
        if feat.strand=='-': children.reverse()
        feat_seq = ''.join(x[0] for x in children)

        print ">%s" % feat.id
        print feat_seq


if __name__ == '__main__':

    from optparse import OptionParser

    p = OptionParser(__doc__)
    p.add_option("--parents", dest="parents", default="mRNA",
            help="list of features to extract, use comma to separate (e.g."
            "'gene,mRNA') [default: %default]")
    p.add_option("--children", dest="children", default="CDS",
            help="list of features to extract, use comma to separate (e.g."
            "'five_prime_UTR,CDS,three_prime_UTR') [default: %default]")

    opts, args = p.parse_args()

    if len(args) != 2:
        sys.exit(p.print_help())
    
    gff_file, fa_file = args
    main(gff_file, fa_file, opts.parents, opts.children)
