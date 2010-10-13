#!/usr/bin/env python
# -*- coding: UTF-8 -*-

r"""
%prog query_fasta subject_fasta [options]

wrapper for calling BLAST in the regions adjacent to given anchors

if qbed and sbed are not specified, then --anchors should point to a file with format:
    qchr\tqstart\tqend\tschr\tsstart\tsend

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
            # it comes with number, line, just want line.
            qbed = dict((k, v[1]) for k, v in qbed.iteritems())
            sbed = dict((k, v[1]) for k, v in sbed.iteritems())

        for line in fh:
            if line[0] == "#": continue
            line = line.rstrip().split("\t")
            q, s = (int(line[1]), int(line[3])) if self.raw else line[:2]
            self.append((qbed[q], sbed[s]))
        self.sort(key=lambda (a, b): (a.seqid, b.seqid, a.start, b.start, a.accn))


    def to_groups(self, distance):
        # not used.
        g = Grouper()
        for name, anchors in itertools.groupby(self, key=lambda a, b: (a.seqid, b.seqid)):
            for ia, qa, sa in enumerate(anchors[:-1]):
                qb, sb = anchors[ia + 1]
                if qb.start - qa.end <= distance and sb.start - sa.end <= distance:
                    g.join((qa, sa), (qb, sb))
        return g

    def _gen_cmds(self, qfasta, sfasta, dist, cmd_template):
        for q, s in self:
            qstart, qend = max(1, q.start - dist - PAD), q.end + dist + PAD
            sstart, send = max(1, s.start - dist - PAD), s.end + dist + PAD
            qtmp = NamedTemporaryFile(delete=False)
            stmp = NamedTemporaryFile(delete=False)
            print >>qtmp, ">%s\n%s" % (q.seqid, str(qfasta[q.seqid][qstart - 1: qend]))
            print >>stmp, ">%s\n%s" % (s.seqid, str(sfasta[s.seqid][sstart - 1: send]))
            qtmp.close()
            stmp.close()
            query_fasta, subject_fasta = qtmp.name, stmp.name
            cmd = cmd_template % locals()
            cmd += "; rm -f %(query_fasta)s %(subject_fasta)s;" % locals()
            # these must match the args accepted by run_blast.
            yield cmd, qstart, sstart, q, s, dist

    def gen_cmds(self, qfasta, sfasta, dist, cmd_template):
        """
        make it easier to parallelize by just yield a list of ncpu commands
        at a time
        """
        n = cpu_count()
        li = []
        for cmd in self._gen_cmds(qfasta, sfasta, dist, cmd_template):
            li.append(cmd)
            if len(li) == n:
                yield li
                li = []
        if li: yield li

class PositionAnchorLine(object):
    __slots__ = ('seqid', 'start', 'end')
    def __init__(self, line):
        self.seqid = line[0]
        self.start = int(line[1])
        self.end   = int(line[2])

    def __str__(self):
        return "%s(%s[%i:%i])" % (self.__class__.__name__, self.seqid,
                           self.start, self.end)
    __repr__ = __str__

class PositionAnchor(Anchor):
    def __init__(self, filename):
        for row in open(filename):
            if row[0] == "#": continue
            row = row.split("\t")
            a = PositionAnchorLine(row[:3])
            b = PositionAnchorLine(row[3:])
            self.append((a, b))
        self.sort(key=lambda (a, b): (a.seqid, b.seqid, a.start, b.start))


def sh(cmd):
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = proc.communicate()
    if err: print >>sys.stderr, err
    return out

def check_locs(locs, bed_line, dist):
    """
    filter out hits that were in the PAD area
    """
    return is_overlapping(locs, bed_line) or \
            calc_dist(locs[0], locs[1], bed_line) <= dist

def is_overlapping(locs, bed_line):
    if locs[1] < locs[0]:
        locs = locs[1], locs[0]
    return locs[1] >= bed_line.start and locs[0] <= bed_line.end


def calc_dist(a, b, bed_line):
    return min(
        abs(a - bed_line.start),
        abs(b - bed_line.start),
        abs(a - bed_line.end),
        abs(b - bed_line.end))

def check_dist(hit, q, s, dist):
    """
    pythagorean yo
    """
    q0, q1 = hit[:2]
    s0, s1 = hit[2:]
    over = is_overlapping((q0, q1), q) and is_overlapping((s0, s1), s)
    if over: return True
    qdist = calc_dist(q0, q1, q)
    sdist = calc_dist(s0, s1, s)
    return ((sdist ** 2 + qdist ** 2) ** 0.5) <= dist

def run_blast(args):
    cmd, qstart, sstart, q, s, dist = args
    data = []
    for line in (l for l in sh(cmd).strip().split("\n") if l):
        if line[0] == "#": continue
        line = line.strip().split("\t")
        # update start, stop positions based on slice.
        line[6:10] = map(int, line[6:10])
        line[6:8] = [l + qstart - 1 for l in line[6:8]]
        if not check_locs(line[6:8], q, dist): continue

        line[8:10] = [l + sstart - 1 for l in line[8:10]]
        # remove those that are outside of the window, but inside the pad.
        if not check_locs(line[8:10], s, dist): continue

        if not check_dist(line[6:10], q, s, dist): continue
        data.append(line)
    return data


def main(qfasta, sfasta, options):
    qfasta = Fasta(qfasta)
    sfasta = Fasta(sfasta)

    if not (options.qbed and options.sbed):
        anchors = PositionAnchor(options.anchors)
    else:
        qbed = Bed(options.qbed)
        sbed = Bed(options.sbed)
        anchors = Anchor(options.anchors, qbed, sbed)
    cpus = cpu_count()
    pool = Pool(cpus)

    for i, command_group in enumerate(anchors.gen_cmds(qfasta, sfasta, options.dist, options.cmd)):
        if not (i - 1) % 500:
            print >>sys.stderr, "complete: %.5f" % (((i - 1.0) * cpus) / len(anchors))
        for lines in pool.map(run_blast, command_group):
            for line in lines:
                line[6:10] = map(str, line[6:10])
                print "\t".join(line)

if __name__ == '__main__':

    import optparse

    p = optparse.OptionParser(__doc__)
    p.add_option("--anchors", dest="anchors",
            help="anchors pair list either as names or integer indexes into"
                      " the bed file")
    p.add_option("--qbed", dest="qbed",
            help="path to qbed")
    p.add_option("--sbed", dest="sbed",
            help="path to sbed")

    p.add_option("--dist", dest="dist",
            default=12000, type="int", help="The extent of flanking regions to search")

    p.add_option("--cmd", dest="cmd", help="the command to be run, but have "
                 "place for 'query_fasta' and 'subject_fasta' as in default "
                  " the command must send the output to STDOUT",
          default="bl2seq -p blastn -e 0.1 -W 15 -D 1 -i %(query_fasta)s -j %(subject_fasta)s")

    (options, fasta_files) = p.parse_args()

    if len(fasta_files) != 2 or not options.anchors:
        sys.exit(p.print_help())
    assert "%(query_fasta)s" in options.cmd
    assert "%(subject_fasta)s" in options.cmd
    qfasta, sfasta = fasta_files

    main(qfasta, sfasta, options)
