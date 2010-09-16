"""
=================
plot simulations
=================

run simulations and plot histograms.

"""

import sys
import os.path as op
from math import floor
import matplotlib; matplotlib.use('Agg')
sys.path.insert(0, op.dirname(__file__))
from fractionation_ga import count_deletion_runs, count_runs, gen_deletions, run_sim
import numpy as np
from matplotlib import pyplot as plt

def plot_runs(runs, pngpath, ax=None):
    """
    given a list of `runs` from count_runs(), plot them
    and save to `pngpath` if given or add to the `ax` if given
    """
    if ax is None:
        f = plt.figure()
        ax = f.add_subplot(1, 1, 1)

    ax.bar([d[0] - 0.5 for d in runs], [d[1] for d in runs], zorder=1)
    ax.set_ylim(0, 20)

    if pngpath is not None:
        ax.set_xlim(xmax=max(d[0] for d in runs))
        ax.set_ylim(ymax=max(d[1] for d in runs))
        plt.savefig(pngpath)

    return ax


def deletion_sim(region_length, deletion_lengths, deletions, simulations=50,
                 count_retentions=False, real_runs=None, plot=True):
    """
    run a simulation of for the given 'region_length' and number of `deletions`
    """
    # initialize as want to know when one is skipped.
    all_runs = dict((x, []) for x in range(400))
    for i in range(simulations):
        ds, runs = gen_deletions(region_length, deletion_lengths, num_deletions=deletions, count_retentions=count_retentions)
        for run_length, n in runs:
            all_runs[run_length].append(n)

    all_runs = dict((k, v) for k, v in all_runs.iteritems() if v != [])

    if plot:
        ranges = {}
        for run_length, counts in sorted(all_runs.items()):
            counts = np.array(counts)
            pct2p5  = int(floor(0.975 * counts.shape[0])) # - 1)
            pct50   = int(floor(0.5 * counts.shape[0])) # - 1)
            pct97p5 = int(floor(0.025 * counts.shape[0])) # - 1)

            if len(counts) != 0:
                counts.sort()
                cmean = counts[pct50]
                c2p5, c97p5 = counts[pct2p5], counts[pct97p5]
                ranges[run_length] = (cmean,  cmean - c2p5, c97p5 - cmean)
            else:
                ranges[run_length] = (0, 0, 0)

        # t, s, e, _, g
        data = np.array([(k, v[0], v[1], v[2]) for k, v in sorted(ranges.items())])
        f = plt.figure()
        a = f.add_subplot(1, 1, 1)
        a.errorbar(data[:, 0], data[:, 1], yerr=np.array([data[:, 2], data[:, 3]]), fmt='+', capsize=300, barsabove=True, zorder=2, ecolor='g')
        a.plot(data[:, 0], data[:, 1], 'ko')
        a.text(0.1, 0.95, str(deletion_lengths), transform=a.transAxes)
        return a, all_runs



def plot(runsdata):
    
    for region, rdict in runsdata.iteritems():
        for homeolog, runs in rdict.items():
            plot_runs(runs, "%s_%s" % (region, homeolog))

def del_sim_from_str(astr, deletion_lengths, simulations=100, count_retentions=False):
    region_length = len(astr)
    nretained = len(astr.replace('_', ''))
    deletions = region_length - nretained

    real_runs = count_deletion_runs(astr)
    #ax = deletion_sim(region_length=region_length,
    return  deletion_sim(region_length=region_length, 
                         deletion_lengths=deletion_lengths,
                         deletions=deletions,
                         simulations=simulations,
                         count_retentions=count_retentions, real_runs=real_runs
                        )
                        

def runall(regions_runs, count_retentions, sim_type, max_del=5):
    """
    run and create a figure for each region, string in d
    """
    deletion_lengths=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    png_template = "data/" + op.splitext(op.basename(regions_runs))[0] \
                + ".%s" % ('retentions' if count_retentions else 'deletions') \
                + ("." + sim_type ) \
                + ((".%i" % max_del) if sim_type == "sim" else "") \
                + ".%s.png"
    for line in open(regions_runs):
        region, runstr = line.strip().split()
        print >>sys.stderr, "saving", png_template % region
        if sim_type == 'one':
            region_with_figure(runstr, png_template % region, count_retentions, deletion_lengths)
        elif sim_type == 'sim':
            sim = run_sim(runstr, max_del)
            region_with_figure(runstr, png_template % region, count_retentions, sim['deletion_lengths'])


def region_with_figure(astr, pngpath, count_retentions, 
                       deletion_lengths):
    """
    plot a figure of the given string
    :param astr: the deletion/retention string
    :param pngpath: where to save the png
    :param count_retentions: bool to count deletions (False) 
                             or retentions (True)
    :param deletion_lengths: list of deletion lenghts for the simulation
                             e.g. [1, 1, 1, 2, 2, 3]
    """
    if count_retentions:
        runs = count_runs(astr)
    else:
        runs = count_deletion_runs(astr)

    actual_fh = open(pngpath.replace('.png','.txt').replace('.sim', '.actual'), 'w')
    print >>sys.stderr, "writing actual to", actual_fh.name
    for run_length, count in runs:
        print >> actual_fh, run_length, count

    ax, runs_from_sim = del_sim_from_str(astr, deletion_lengths, count_retentions=count_retentions)
    txt_fh = open(pngpath.replace('.png', '.txt'), 'w')
    print >>sys.stderr, "writing sims to", txt_fh.name
    for run_length in runs_from_sim:
        print >>txt_fh, run_length, sum(runs_from_sim[run_length])
    plot_runs(runs, None, ax=ax)
    ax.set_xlim(xmax=max(x[0] for x in runs))
    ax.set_ylim(ymax=1.5 * max(x[1] for x in runs))
    ax.set_xlabel('run length')
    ax.set_ylabel('count')
    if count_retentions:
        ax.set_title('retentions')
    else:
        ax.set_title('deletions')
    plt.savefig(pngpath)
    

if __name__ == "__main__":

    if sys.argv[1].endswith('.txt'):
        runs_file = sys.argv[1]
        do_count_retentions = int(sys.argv[2])
        sim_type = sys.argv[3]
        max_del = len(sys.argv) > 3 and int(sys.argv[4]) or 1
        assert sim_type in ('sim', 'one', 'gtest')
        runall(runs_file, do_count_retentions, sim_type=sim_type, max_del=max_del)

    """
        sys.exit()
    # SEE region_with_figure for simpler entry point.
    astr = sys.argv[1] 
    # send in name for plot.
    name = sys.argv[2] if len(sys.argv) > 2 else None  
    count_retentions = len(sys.argv) > 3

    deletion_lengths=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
    ax = del_sim_from_str(astr, deletion_lengths, count_retentions=count_retentions)
    if count_retentions:
        runs = count_runs(astr)
    else:
        runs = count_deletion_runs(astr)

    print runs
    plot_runs(runs, name, ax=ax)
    ax.set_xlim(xmax=max(x[0] for x in runs))
    ax.set_ylim(ymax=max(x[1] for x in runs))
    ax.set_xlim(xmax=50)
    ax.set_ylim(ymax=110)
    ax.set_xlabel('run length')
    ax.set_ylabel('count')
    ax.set_title(name)
    plt.savefig('/var/www/t/sb_zm_fractionation/2497/bar_%s.png' % (name, ))
    """
