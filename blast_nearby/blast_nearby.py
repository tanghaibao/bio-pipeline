#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog query_fasta subject_fasta [options]

wrapper for calling BLAST in the regions adjacent to given anchors

"""


import sys
from pyfasta import Fasta
from bed_utils import Bed
from grouper import Grouper
import itertools
from tempfile import NamedTemporaryFile
from multiprocessing import Pool, cpu_count
from subprocess import Popen, PIPE

PAD = 1000

class Anchor(list):
    def __init__(self, filename, qbed, sbed):
        self.filename = filename
        fh = open(self.filename)
        while True:
            line = fh.readline().strip()
            if line[0] == "#": continue
            line = line.split("\t")
            # raw format
            self.raw = len(line) > 2 and line[1].isdigit() and line[3].isdigit()
            fh.seek(0)
            break

        # index by name.
        if not self.raw:
            qbed, sbed = qbed.get_order(), sbed.get_order()

        for line in fh:
            if line[0] == "#": continue
            line = line.split("\t")
            q, s = (int(line[1]), int(line[3])) if self.raw else line[:2]
            self.append((qbed[q], sbed[s]))
        self.sort(key=lambda (a, b): (a.seqid, b.seqid, a.start, b.start, a.accn))


    def to_groups(self, distance):
        g = Grouper()
        for name, anchors in itertools.groupby(self, key=lambda a, b: (a.seqid, b.seqid)):
            for ia, qa, sa in enumerate(anchors[:-1]):
                qb, sb = anchors[ia + 1]
                if qb.start - qa.end <= distance and sb.start - sa.end <= distance:
                    g.join((qa, sa), (qb, sb))
        return g

    def _gen_cmds(self, qfasta, sfasta, dist, params):
        for q, s in self:
            qstart, qend = max(1, q.start - dist - PAD), q.end + dist + PAD
            sstart, send = max(1, s.start - dist - PAD), s.end + dist + PAD
            qtmp = NamedTemporaryFile(delete=False)
            stmp = NamedTemporaryFile(delete=False)
            print >>qtmp, ">%s\n%s" % (q.seqid, str(qfasta[q.seqid][qstart - 1: qend]))
            print >>stmp, ">%s\n%s" % (s.seqid, str(qfasta[q.seqid][qstart - 1: qend]))
            qtmp.close()
            stmp.close()
            cmd = "bl2seq -p blastn -W 15 -D 1 -i %s -j %s; rm -f %s %s"
            cmd %= (qtmp.name, stmp.name, qtmp.name, stmp.name)
            yield cmd, qstart, sstart

    def gen_cmds(self, qfasta, sfasta, dist, params):
        n = cpu_count()
        li = []
        for cmd in self._gen_cmds(qfasta, sfasta, dist, params):
            li.append(cmd)
            if len(li) == n:
                yield li
                li = []
        if li: yield li

def sh(cmd):
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = proc.communicate()
    if err: print >>sys.stderr, err
    return out

def run_blast(args):
    #print >>sys.stderr, args
    cmd, qstart, sstart = args
    data = []
    for line in sh(cmd).strip().split("\n"):
        #print >>sys.stderr, line
        if line[0] == "#": continue
        line = line.strip().split("\t")
        # update start, stop positions based on slice.
        line[6:10] = map(int, line[6:10])
        line[6:8] = [l + qstart - 1 for l in line[6:8]]
        line[8:10] = [l + sstart - 1 for l in line[8:10]]
        line[6:10] = map(str, line[6:10])
        data.append("\t".join(line))
    return data


def main(qfasta, sfasta, options):
    qfasta = Fasta(qfasta)
    sfasta = Fasta(sfasta)

    qbed = Bed(options.qbed)
    sbed = Bed(options.qbed)

    anchors = Anchor(options.anchors, qbed, sbed)
    # TODO: use shlex to parse options.params. then dict.update
    pool = Pool(cpu_count())
    for i, command_group in enumerate(anchors.gen_cmds(qfasta, sfasta, options.dist, options.parameters)):
        if not (i - 1) % 1000:
            print >>sys.stderr, "complete: %.5f" % (((i - 1) / 8) / len(anchors))
        for lines in pool.map(run_blast, command_group):
            print "\n".join(lines)

if __name__ == '__main__':

    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--anchors", dest="anchors",
            help="anchors pair list either as names or integer indexes into"
                      " the bed file")
    parser.add_option("--qbed", dest="qbed",
            help="path to qbed")
    parser.add_option("--sbed", dest="sbed",
            help="path to sbed")

    parser.add_option("--dist", dest="dist",
            default=12000, type="int", help="The extent of flanking regions to search")
    parser.add_option("--parameters", dest="parameters",
            default=" -m 8 -p blastn ", help="parameters to run.")

    (options, fasta_files) = parser.parse_args()

    if len(fasta_files) != 2 or not \
            all((options.qbed, options.sbed, options.anchors)):
        sys.exit(parser.print_help())
    qfasta, sfasta = fasta_files

    main(qfasta, sfasta, options)

