#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog folder_of_mafs > out.bed

convert a folder of maf alignments to the bed features 
then useful to check coverage, etc.
"""

import os
import sys
import numpy as np

from glob import glob
from bx.align import maf


def hamming_distance(s1, s2):
    """
    >>> hamming_distance("attttgggaa", "ATTGagaGAT")
    4
    """
    assert len(s1)==len(s2), ",".join((s1, s2))
    s1 = np.array(s1.upper(), dtype="c")
    s2 = np.array(s2.upper(), dtype="c")
    return np.sum(s1 != s2)


def main(folder, options):
    flist = glob(folder + "*.maf")
    print >>sys.stderr, folder, "contains %d files" % len(flist) 
    j = 0
    for f in flist:
        try:
            reader = maf.Reader(file(f))
        except:
            print >>sys.stderr, "[warning] %s not readable" % f
            continue
        for rec in reader:
            a, b = rec.components
            length = len(a.text) 
            distance = hamming_distance(a.text, b.text) * 100. / length

            for a, atag in zip((a, b), "ab"):
                print "\t".join(str(x) for x in (a.src, a.forward_strand_start, \
                        a.forward_strand_end, \
                        "%s%07d%s" % (folder, j, atag), \
                        length, distance, \
                        a.text))
                sys.stdout.flush()

            j += 1


if __name__ == '__main__':

    import doctest
    doctest.testmod()

    from optparse import OptionParser
    parser = OptionParser(usage=__doc__)

    (options, args) = parser.parse_args()

    if len(args)==1: 
        folder = args[0]
    else:
        sys.exit(parser.print_help())

    main(folder, options)
