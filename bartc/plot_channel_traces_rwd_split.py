"""
Given a tuple of subject, dataset, line plot peri-event power in each
LFP channel
Split this up based on reward.
Based on plot_channel_traces.py
"""
from __future__ import division
import numpy as np
import pandas as pd
import physutils
from hephys import dbio
import os
import itertools

from plot_channel_traces import *

def get_traces_split(dtup, event, bands):

    # load data
    os.chdir(os.path.expanduser('~/code/hephys/bartc'))
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

    # get lfp data
    print "Fetching Data..."
    lfp = dbio.fetch_all_such_LFP(dbname, *dtup)
        
    # bandpass filter
    print "Filtering..."
    lfp = lfp.bandlimit(bands)

    # decimate to 100 Hz effective sampling
    print "Decimating..."
    lfp = lfp.decimate(5)

    # instantaneous power
    print "Calculating Power..."
    lfp = lfp.instpwr()

    # remove censored regions
    print "Censoring..."
    lfp = lfp.censor()

    # get events
    evt = dbio.fetch(dbname, 'events', *dtup)
    evtseries = evt[event].dropna()

    # median split by inflate time
    infl_times = evt['inflate_time']
    med_inflate_time = np.median(infl_times)
    evt_grps = []
    evt_grps.append(evtseries[infl_times <= med_inflate_time])
    evt_grps.append(evtseries[infl_times > med_inflate_time])
    grpnames = ['Low value', 'High Value']

    subsets = []
    for grp in evt_grps:
        # split lfp around stops
        lfp_split = lfp.evtsplit(grp, Tpre, Tpost)

        # group by time and get median for each channel for each time
        medians = lfp_split.groupby(level=1).median()

        # get median across channels
        grand_median = medians.median(axis=1)

        subsets.append(grand_median)

    combined = pd.concat(subsets, axis=1)
    combined.columns = grpnames

    # make peri-stop frame
    df = physutils.LFPset(combined, meta=lfp.meta.copy()).zscore()
    df = df.smooth(smwid)

    return df

if __name__ == '__main__':
    evtnames = ['banked', 'popped']
    bandlist = [['theta'], ['alpha'], ['beta']]
    tuplist = [(17, 2), (18, 1), (20, 1), (22,1), (23, 1), (30, 1)]
    for (dtup, event, bands) in itertools.product(tuplist, evtnames, bandlist):
        print dtup, event, bands
        df = get_traces_split(dtup, event, bands)
        print_from_R(df, event, bands)