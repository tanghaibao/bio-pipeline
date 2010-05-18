#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog query_fasta subject_fasta [options]

wrapper for calling BLAST in the regions adjacent to given anchors

"""


import sys

def main(fasta_files, options):
    pass

if __name__ == '__main__':
    
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--anchors", dest="anchors",
            help="anchors pair list")
    parser.add_option("--qbed", dest="qbed", 
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed", 
            help="path to sbed")

    parser.add_option("--dist", dest="dist",
            default=12000, type="int", help="The extent of flanking regions to search")

    (option, fasta_files) = parser.parse_args()

    if not len(fasta_files)==2:
        sys.exit(parser.print_help())

    main(fasta_files, options)

