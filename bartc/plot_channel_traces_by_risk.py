
# # Plot average peri-event power for all risk conditions

from __future__ import division
import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas as pd
import physutils
from hephys import dbio
import os
import itertools

# set up relevant params
Tpre = -2  # time before event to grab
Tpost = 0.5  # time after event to grab
smwid = 0.2  # smoothing window (in s)

def get_traces(dtup, event, bands):

    # load data
    os.chdir(os.path.expanduser('~/code/hephys/bartc'))
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')
    lfp = dbio.fetch_all_such_LFP(dbname, *dtup)


    # bandpass filter
    lfp = lfp.bandlimit(bands)

    # decimate to 100 Hz
    lfp = lfp.decimate(5)

    # instantaneous power
    lfp = lfp.instpwr()

    # censor
    lfp = lfp.censor()

    # zscore to facilitate cross-channel comparisons
    lfp = lfp.zscore()

    # get events
    evt = dbio.fetch(dbname, 'events', *dtup)

    # restric to non-control trials
    evt = evt[evt.trial_type.isin([1, 2, 3])]

    # get only those trials where outcome is not NA
    evtseries = evt[[event, 'trial_type']].dropna()

    # split lfp around stops
    lfp_split = lfp.evtsplit(evtseries[event], Tpre, Tpost)

    # groupby accepts a list of functions to group by
    # each function gets an index tuple
    group_by_trial_type = lambda x: evtseries.iloc[x[0]].trial_type
    group_by_time = lambda x: x[1]
    grpfuns = [group_by_trial_type, group_by_time]
    grouped = lfp_split.groupby(grpfuns)

    # now get median across all trials for each (type, time, channel)
    med_by_type = grouped.median()
    med_across_chans = med_by_type.median(axis=1)

    # make trial type a column
    medians = med_across_chans.unstack(level=0)
    medians.index.name = 'time'
    medians.columns = map(int, medians.columns)

    # make peri-stop frame
    df = physutils.LFPset(medians, meta=lfp.meta.copy()).zscore()
    df = df.smooth(smwid)

    return df


def make_filename(df, name, bands):
    base = '~/Dropbox/hephys/media/figs/' 
    pieces = list(df.meta['tuple']) + ['_'.join(bands), name, 'chanplot', 'pdf']
    return "\'" + base + '.'.join(map(str, pieces)) + "\'"

def print_from_R(df, name, bands):
    # prepare dataframe for passing to R
    rdat = df.copy()
    rdat = rdat.reset_index()
    rdat = pd.melt(rdat, id_vars='time')
    rdf = com.convert_to_r_dataframe(rdat)

    fname = make_filename(df, name, bands)

    # load up R
    R = robjects.r
    R("""source('helpers.R')""")

    # make plot
    R("pdf(file=" + fname + ", paper='USr', width=11, height=8.5)")

    pp = R['chanmeans'](rdf)
    R.plot(pp)

    R('dev.off()')

if __name__ == '__main__':
    evtnames = ['banked', 'popped']
    bandlist = [['theta'], ['alpha'], ['beta']]
    tuplist = [(17, 2), (18, 1), (20, 1), (22,1), (23, 1), (30, 1)]
    for (dtup, event, bands) in itertools.product(tuplist, evtnames, bandlist):
        print dtup, event, bands
        df = get_traces(dtup, event, bands)
        print_from_R(df, event, bands)

