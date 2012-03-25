#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog -i query.fa -d database.fa [options]

run LASTZ similar to the BLAST interface, and generates -m8 tabular format
"""

import os
import sys
import math
import logging
logging.basicConfig(level=logging.DEBUG)

from subprocess import Popen, PIPE
from multiprocessing import Process, cpu_count, Lock


blast_fields = "query,subject,pctid,hitlen,nmismatch,ngaps,"\
        "qstart,qstop,sstart,sstop,evalue,score"

lastz_fields = "name2,name1,identity,nmismatch,ngap,"\
        "start2+,end2+,strand2,start1,end1,strand1,score"

# conversion between blastz and ncbi is taken from Kent src
# src/lib/blastOut.c
# this is not rigorous definition of e-value (assumes human genome) !!
blastz_score_to_ncbi_bits = lambda bz_score: bz_score * 0.0205

def blastz_score_to_ncbi_expectation(bz_score):
    bits = blastz_score_to_ncbi_bits(bz_score)
    log_prob = -bits * 0.693147181
    # this number looks like.. human genome?
    return 3.0e9 * math.exp(log_prob)


def lastz_to_blast(row):
    # conver the lastz tabular to the blast tabular, see headers above
    atoms = row.strip().split("\t")
    name1, name2, coverage, identity, nmismatch, ngap, \
            start1, end1, strand1, start2, end2, strand2, score = atoms
    identity = identity.replace("%", "")
    hitlen = coverage.split("/")[1]
    score = float(score)
    same_strand = (strand1 == strand2)
    if not same_strand:
        start2, end2 = end2, start2

    evalue = blastz_score_to_ncbi_expectation(score)
    score = blastz_score_to_ncbi_bits(score)
    evalue, score = "%.2g" % evalue, "%.1f" % score
    return "\t".join((name1, name2, identity, hitlen, nmismatch, ngap, \
            start1, end1, start2, end2, evalue, score))


def lastz(k, n, bfasta_fn, out_fh, lock, lastz_path, extra):
    lastz_bin = lastz_path or "lastz"

    lastz_cmd = "%s --format=general-:%s "\
            "--ambiguous=iupac %s[multiple,unmask,nameparse=darkspace]"\
            " %s[unmask,nameparse=darkspace,subsample=%d/%d] %s"
    lastz_cmd %= (lastz_bin, lastz_fields, bfasta_fn, afasta_fn, k, n, extra)

    proc = Popen(lastz_cmd, bufsize=1, stdout=PIPE, shell=True)

    logging.debug("job <%d> started: %s" % (proc.pid, lastz_cmd))
    for row in proc.stdout:
        brow = lastz_to_blast(row)
        lock.acquire()
        print >>out_fh, brow
        out_fh.flush()
        lock.release()
    logging.debug("job <%d> finished" % proc.pid)


def main(options, afasta_fn, bfasta_fn, out_fh, extra):

    lastz_path = options.lastz_path
    # split on query so check query fasta sequence number
    afasta_num = sum(1 for x in open(afasta_fn) if x[0]=='>')
    cpus = min(options.cpus, cpu_count(), afasta_num)
    logging.debug("Dispatch job to %d cpus" % cpus)

    processes = []
    lock = Lock()
    for k in xrange(cpus):
        pi = Process(target=lastz, args=(k+1, cpus, bfasta_fn, out_fh, lock,
            lastz_path, extra))
        pi.start()
        processes.append(pi)

    for pi in processes:
        pi.join()


if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser(__doc__)
    parser.add_option("-i", dest="query",
            help="query sequence file in FASTA format")
    parser.add_option("-d", dest="target",
            help="database sequence file in FASTA format")
    parser.add_option("-o", dest="outfile",
            help="BLAST output [default: stdout]")
    parser.add_option("-a", "-A", dest="cpus", default=1, type="int",
            help="parallelize job to multiple cpus [default: %default]")
    parser.add_option("--path", dest="lastz_path", default=None,
            help="specify LASTZ path")
    parser.add_option("--lastz-params", dest="extra", default="",
            help="pass in LASTZ parameter string (please quote the string)")

    (options, args) = parser.parse_args()

    try:
        afasta_fn = options.query
        assert os.path.exists(afasta_fn), ("%s does not exist" % afasta_fn)
        bfasta_fn = options.target
        assert os.path.exists(bfasta_fn), ("%s does not exist" % bfasta_fn)
        out_fh = file(options.outfile, "w") if options.outfile else sys.stdout
    except Exception, e:
        print >>sys.stderr, str(e)
        sys.exit(parser.print_help())

    if not all((afasta_fn, bfasta_fn)):
        sys.exit(parser.print_help())

    main(options, afasta_fn, bfasta_fn, out_fh, options.extra)

