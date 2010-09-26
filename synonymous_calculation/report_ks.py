#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
%prog <ks file>

generate a report given a Ks result file (as produced by synonymous_calc.py).
describe the median Ks, Ka values, as well as the distribution in stem-leaf plot
'''

import sys
import csv
import numpy as np
from collections import namedtuple

fields = "pair yn_ks yn_ka ng_ks ng_ka"
descriptions = {
        'pair': 'Gene pair',
        'yn_ks': 'Yang-Nielson method of Ks estimate',
        'yn_ka': 'Yang-Nielson method of Ka estimate',
        'ng_ks': 'Nei-Gojobori method of Ks estimate',
        'ng_ka': 'Nei-Gojobori method of Ka estimate'}

KsLine = namedtuple("KsLine", fields)

def read_ks_file(ks_file):
    reader = csv.reader(open(ks_file, "rb"))
    reader.next() # header
    data = []
    for row in reader:
        for i, a in enumerate(row):
            if i==0: continue
            row[i] = float(row[i])
        data.append(KsLine._make(row))
    return data 


def stem_leaf_plot(data, bins=np.arange(0, 1.1, .1)):
    '''
    generate stem and leaf plot given a collection of numbers
    '''
    hist, bin_edges = np.histogram(data, bins=bins)
    width = 50 # the width (height) of the distribution
    hist = hist * width / hist.sum()
    for b, h in zip(bin_edges, hist):
        print "%s|%s" % (("%.1f" % b).rjust(5), "=" * h)


def main(ks_file):

    data = read_ks_file(ks_file)
    print 'File `%s` contains a total of %d gene pairs' % (ks_file, len(data))
    print '-' * 60
    for f in fields.split()[1:]:
        columndata = [getattr(x, f) for x in data]
        print descriptions[f], ':', np.median(columndata)
        bins = np.arange(0, 2.2, .2) if 'ks' in f else np.arange(0, .6, .1)
        stem_leaf_plot(columndata, bins=bins)


if __name__ == '__main__':
    
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==1:
        ks_file = args[0]
    else:
        sys.exit(p.print_help())

    main(ks_file)
