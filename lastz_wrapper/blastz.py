#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog -i query.fa -d database.fa [options]

run LASTZ similar to the BLAST interface, and generates -m8 tabular format
"""

import os
import sys
import math
import shutil
import tempfile
import logging
logging.basicConfig(level=logging.DEBUG)

from subprocess import Popen, PIPE
from multiprocessing import Process, cpu_count, Lock
from Bio import SeqIO
from random import shuffle


blast_fields = "query,subject,pctid,hitlen,nmismatch,ngaps,"\
        "qstart,qstop,sstart,sstop,evalue,score"

lastz_fields = "name2,name1,identity,nmismatch,ngap,"\
        "start2,end2,start1,end1,score"

# conversion between blastz and ncbi is taken from Kent src
# src/lib/blastOut.c
# this is not rigorous definition of e-value (assumes human genome) !!
blastz_score_to_ncbi_bits = lambda bz_score: bz_score * 0.0205

def blastz_score_to_ncbi_expectation(bz_score):
    bits = blastz_score_to_ncbi_bits(bz_score)
    log_prob = -bits * 0.301029996
    # this number looks like.. human genome?
    return 3.0e9 * math.exp(log_prob)


def lastz_to_blast(row):
    # conver the lastz tabular to the blast tabular, see headers above
    atoms = row.strip().split("\t")
    name1, name2, coverage, identity, nmismatch, ngap, \
            start1, end1, start2, end2, score = atoms
    identity = identity.replace("%", "")
    hitlen = coverage.split("/")[1]
    score = float(score)

    evalue = blastz_score_to_ncbi_expectation(score)
    score = blastz_score_to_ncbi_bits(score)
    evalue, score = "%.2g" % evalue, "%.1f" % score
    return "\t".join((name1, name2, identity, hitlen, nmismatch, ngap, \
            start1, end1, start2, end2, evalue, score))


def lastz(afasta_fn, bfasta_fn, out_fh, lock, lastz_path):
    lastz_bin = "lastz" if lastz_path is None else lastz_path

    lastz_cmd = "%s --format=general-:%s --ambiguous=iupac %s[multiple,unmask,nameparse=darkspace] %s[unmask,nameparse=darkspace]"
    lastz_cmd %= (lastz_bin, lastz_fields, bfasta_fn, afasta_fn)
    logging.debug(lastz_cmd)

    proc = Popen(lastz_cmd, bufsize=1, stdout=PIPE, shell=True)

    logging.debug("job <%d> started" % proc.pid)
    for row in proc.stdout:
        brow = lastz_to_blast(row)
        lock.acquire()
        print >>out_fh, brow
        out_fh.flush()
        lock.release()
    logging.debug("job <%d> finished" % proc.pid)


def chunks(L, n):
    # yield successive n-sized chunks from list L
    # shuffle since largest chrs are often first, get bad
    # distribution to split files.
    shuffle(L)
    for i in xrange(0, len(L), n):
        yield L[i:i+n]


def newnames(oldname, n):
    return ["%s.%02d" % (oldname, i) for i in xrange(n)]


def main(options, afasta_fn, bfasta_fn, out_fh):

    lastz_path = options.lastz_path
    cpus = min(options.cpus, cpu_count())
    logging.debug("Dispatch job to %d cpus" % cpus)

    if cpus > 1:
        # split input fasta into chunks
        print >>sys.stderr, "splitting fasta to %i files" % cpus
        recs = list(SeqIO.parse(afasta_fn, "fasta"))
        temp_dir = tempfile.mkdtemp(prefix=".split", dir=os.getcwd())
        names = newnames(os.path.join(temp_dir, os.path.basename(afasta_fn)), cpus)
        logging.debug("Temporary directory %s created for %s" % \
                (temp_dir, afasta_fn))
        # TODO: make sure this works for all cases
        # zip only exhausts shorted
        chunk_size = len(recs) / cpus + 1
        nchunks = 0
        for name, chunk in zip(names, chunks(recs, chunk_size)):
            if chunk is None: break
            logging.debug("NAME:" + name)
            nchunks += len(chunk)
            SeqIO.write(chunk, name, "fasta")
        assert nchunks == len(recs), (nchunks, len(recs))
    else:
        names = [afasta_fn]

    processes = []
    lock = Lock()
    for name in names:
        if not os.path.exists(name): continue
        pi = Process(target=lastz, args=(name, bfasta_fn, out_fh, lock, lastz_path))
        pi.start()
        processes.append(pi)

    for pi in processes:
        pi.join()

    if cpus > 1:
        shutil.rmtree(temp_dir)


if __name__ == '__main__':

    from optparse import OptionParser

    parser = OptionParser(__doc__)
    parser.add_option("-i", dest="query",
            help="query sequence file in FASTA format")
    parser.add_option("-d", dest="target",
            help="database sequence file in FASTA format")
    parser.add_option("-o", dest="outfile",
            help="BLAST output [default: stdout]")
    parser.add_option("-A", dest="cpus", default=1, type="int",
            help="parallelize job to multiple cpus [default: %default]")
    parser.add_option("--path", dest="lastz_path", default=None,
            help="specify LASTZ path")

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

    main(options, afasta_fn, bfasta_fn, out_fh)

