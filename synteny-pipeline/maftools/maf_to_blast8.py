#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
%prog maf_folder > out.blast

From a folder of .maf files, generate a single .blast file that is the BLAST m8
format
'''

import sys
from glob import glob

from bx.align import maf
from blastz import blastz_score_to_ncbi_expectation, blastz_score_to_ncbi_bits

# same as in blastz.py
blast_fields = "query,subject,pctid,hitlen,nmismatch,ngaps,"\
                "qstart,qstop,sstart,sstop,evalue,score"


def alignment_details(a, b):
    # TODO: not implemented
    assert len(a)==len(b), "seq len do not match"
    #a = np.array(a.upper(), dtype='c')
    #b = np.array(b.upper(), dtype='c')
    #matches = np.sum(a == b)
    pctid = 0
    nmismatch = 0
    ngaps = 0
    return pctid, nmismatch, ngaps


def maf_to_blast8(f):
    try:
        print >>sys.stderr, "reading %s" % f
        reader = maf.Reader(open(f))
    except:
        print >>sys.stderr, "[warning] %s not readable" % f
        return
    
    for rec in reader:
        a, b = rec.components
        query = a.src
        subject = b.src
        qstart = a.forward_strand_start
        qstop = a.forward_strand_end
        sstart = b.forward_strand_start
        sstop = b.forward_strand_end
        score = rec.score

        evalue = blastz_score_to_ncbi_expectation(score)
        score = blastz_score_to_ncbi_bits(score)
        evalue, score = "%.2g" % evalue, "%.1f" % score
        hitlen = len(a.text)

        #print a.text
        #print b.text
        pctid, nmismatch, ngaps = alignment_details(a.text, b.text)
        print "\t".join(str(x) for x in (query,subject,pctid,hitlen,nmismatch,ngaps,
                qstart,qstop,sstart,sstop,evalue,score))


def main(maf_folder):
    flist = sorted(glob(maf_folder + "*.maf"))
    print >>sys.stderr, maf_folder, "contains %d files" % len(flist)

    for f in flist:
        maf_to_blast8(f)


if __name__ == '__main__':

    from optparse import OptionParser

    p = OptionParser(__doc__)
    (opts, args) = p.parse_args()

    try:
        maf_folder = args[0]
    except Exception, e:
        print >>sys.stderr, str(e)
        sys.exit(p.print_help())

    main(maf_folder)
