#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
from optparse import OptionParser

from bx.align import maf

def main(options, args):
    infile, chr1, chr2 = args

    in_file = args[0]
    base, ext = os.path.splitext(in_file)
    out_file = "%(base)s.%(chr1)s_vs_%(chr2)s_filtered%(ext)s" % locals()

    fp = file(in_file)
    fw = file(out_file, "w")

    reader = maf.Reader(fp)
    writer = maf.Writer(fw)
    for rec in reader:
        c1, c2 = rec.components[0].src, rec.components[1].src
        if (chr1, chr2) == (c1, c2) or (chr1, chr2) == (c2, c1):
            writer.write(rec)

if __name__ == '__main__':
    usage = "usage: %prog [options] infile chr1 chr2"
    parser = OptionParser(usage=usage, version="%prog 1.0")
    
    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit()

    main(options, args)

