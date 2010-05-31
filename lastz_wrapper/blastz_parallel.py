#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog -i query.fa -d database.fa [options]

run LASTZ similar to the BLAST interface, and generates -m8 tabular format
"""

import os
import sys
import shutil
import tempfile
import logging
logging.basicConfig(level=logging.ERROR)

from subprocess import Popen, PIPE
from multiprocessing import Process, cpu_count
from Bio import SeqIO


blast_fields = "query,subject,pctid,hitlen,nmismatch,ngaps,"\
        "qstart,qstop,sstart,sstop,evalue,score"

lastz_fields = "name2,name1,identity,nmismatch,ngap,"\
        "start2,end2,start1,end1,score"


def lastz_to_blast(row):
    # conver the lastz tabular to the blast tabular, see headers above
    atoms = row.strip().split("\t")
    name1, name2, coverage, identity, nmismatch, ngap, \
            start1, end1, start2, end2, score = atoms
    identity = identity.replace("%", "")
    hitlen = coverage.split("/")[1]
    # TODO: evalue formula lookup
    evalue = str(1e-50)
    return "\t".join((name1, name2, identity, hitlen, nmismatch, ngap, \
            start1, end1, start2, end2, evalue, score))


def lastz(afasta_fn, bfasta_fn, out_fh):
    lastz_cmd = "lastz --format=general-:%s --ambiguous=iupac %s[multiple] %s"
    lastz_cmd %= (lastz_fields, bfasta_fn, afasta_fn) 
    logging.debug(lastz_cmd)

    proc = Popen(lastz_cmd, bufsize=1, stdout=PIPE, stderr=PIPE, shell=True)

    logging.debug("job <%d> started" % proc.pid)
    for row in proc.stdout:
        print >>out_fh, lastz_to_blast(row)
    logging.debug("job <%d> finished" % proc.pid)


def lastz_process(afasta_fn, bfasta_fn, out_fh):
    lastz(afasta_fn, bfasta_fn, out_fh)


def chunks(L, n):
    # yield successive n-sized chunks from list L 
    for i in xrange(0, len(L), n):
        yield L[i:i+n]


def newnames(oldname, n):
    return ["%s.%02d" % (oldname, i) for i in xrange(n)]


def main(cpus, afasta_fn, bfasta_fn, out_fh):
    cpus = min(cpus, cpu_count())
    logging.debug("Dispatch job to %d cpus" % cpus)

    # split input fasta into chunks
    recs = list(SeqIO.parse(afasta_fn, "fasta"))
    temp_dir = tempfile.mkdtemp(prefix=".split", dir=os.getcwd())
    names = newnames(os.path.join(temp_dir, afasta_fn), cpus)
    logging.debug("Temporary directory %s created for %s" % \
            (temp_dir, afasta_fn))
            
    chunk_size = len(recs) / cpus + 1
    for name, chunk in zip(names, chunks(recs, chunk_size)):
        if chunk is None: break 
        SeqIO.write(chunk, name, "fasta")

    processes = []
    for name in names:
        pi = Process(target=lastz_process, 
                     args=(name, bfasta_fn, out_fh))
        pi.start()
        processes.append(pi)

    for pi in processes:
        pi.join()

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
            help="serialize job to multiple cpus [default: 1]")

    (options, args) = parser.parse_args()

    try:
        afasta_fn = options.query
        bfasta_fn = options.target
        out_fh = file(options.outfile, "w") if options.outfile else sys.stdout
    except:
        sys.exit(parser.print_help())

    if not all((options.query, options.target)):
        sys.exit(parser.print_help())

    main(options.cpus, afasta_fn, bfasta_fn, out_fh)

