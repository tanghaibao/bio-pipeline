#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
[usage] %prog block_name
To test for nested models in the simulation of gaps after polyploidy:

estimate the likelihood for model 1, 2, 3, .. (all the files need to be in
`sims/` folder)
- actual observations in 'actual.1.A.txt'
- model 1 (allowing deletion size of 1) in 'sim.1.A.txt'
- model 2 (allowing deletion size of 1, 2) in 'sim.2.A.txt'

each model will allow one more deletion size in brentp's genetic algorithm
see <http://code.google.com/p/bpbio/source/browse/trunk/scripts/fractionation/>

Likelihood ratio test follows to see which model gives the best likelihood while
keeping the model as simple as possible (Occam's Razer) '''

import os
import sys
import numpy as np

from scipy.stats import chisqprob

prefix = "runs.a.deletions."
N = 10 # ignore deletion size >N 


class Model(object):

    def __init__(self, name, counts, probs, df):
        self.name = name
        self.counts = counts
        self.probs = probs
        self.df = df
        self.lnL = self.log_likelihood(counts, probs)

    def __str__(self):
        return '%s\n%s\nlnL = %.1f' % (self.name, self.probs, self.lnL)

    def log_likelihood(self, counts, probs):
        return np.add.reduce(counts * np.log(probs))


def read_distribution(filename, freq=False, N=N):
    '''
    file look like this:
    1 100
    2 30
    ...
    '''
    fp = open(filename)
    data = [int(row.split()[1]) for row in fp]
    a = np.array(data[:N], 'f')
    if freq:
        a /= a.sum()
    return a


def likelihood_ratio_test(counts, model1, model2):
    # see formula <http://en.wikipedia.org/wiki/Likelihood-ratio_test>
    print 'Test %s and %s' % (model1.name, model2.name)
    D = -2 * (model1.lnL - model2.lnL)
    df = model2.df - model1.df
    p_value = chisqprob(D, df) 

    print 'D = %.1f, df = %d, P-value = %.2g' % (D, df, p_value)


def print_banner(char='-', N=80):
    print char * N


def pairwise(iterable):
    '''
    itertools recipe:
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    '''
    from itertools import tee, izip
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def main(block='A'):

    os.chdir('sims')

    print 'Summarizing models for block %s' % block
    print_banner('=')
    counts = read_distribution(prefix + 'actual.1.%s.txt' % block)
    models = []
    for i in range(1, 5):
        probs = read_distribution(prefix + 'sim.%d.%s.txt' % (i, block), freq=True)
        model = Model('Model%d' % i, counts, probs, i)
        models.append(model)

    print 'Observed counts'
    print counts 
    print_banner()
    for model in models:
        print model
        print_banner()

    for model1, model2 in pairwise(models):
        likelihood_ratio_test(counts, model1, model2)
        print_banner()


if __name__ == '__main__':

    from optparse import OptionParser

    p = OptionParser(__doc__)
    if len(sys.argv) != 2:
        sys.exit(p.print_help())

    main(block=sys.argv[1])
