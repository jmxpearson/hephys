import numpy as np
import pandas as pd
from physutils import *
from physclasses import *

# define some useful numbers
np.random.seed(12345)

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/lfplist')

# group by (patient, dataset) entries:
groups = setlist.groupby(['patient', 'dataset'])

# iterate over groups
for name, grp in groups:

    allchans = []

    # iterate over channels within groups
    for ind, series in grp.iterrows():

        dtup = tuple(series.values)

        print dtup

        # read in data
        print 'Reading LFP...'
        lfp = fetch_LFP(dbname, *dtup)

        # de-mean
        print 'Removing mean...'
        lfp = lfp.demean()

        # break out by frequency bands
        print 'Filtering by frequency...'
        filters = ['delta', 'theta', 'alpha']
        banded = lfp.bandlimit(filters)

        # decimate down to 40 Hz
        print 'Decimating...'
        banded = banded.decimate(5)

        # get instantaneous power
        print 'Calculating power...'
        banded = banded.instpwr()

        # handle censoring
        print 'Censoring...'
        banded = banded.censor()

        # standardize per channel
        print 'Standardizing regressors...'
        banded = banded.zscore()


        # append to whole dataset
        allchans.append(banded)

    # concatenate data from all channels
    print 'Merging channels...'
    groupdata = pd.concat(allchans)

    # specify peri-event times
    dt = 1. / np.array(banded.meta['sr']).round(3)  # dt in ms
    Tpre = 2  # time relative to event to start
    Tpost = 1.5  # time following event to exclude

    # running average of power
    print 'Running mean...'
    meanpwr = pd.rolling_mean(groupdata, np.ceil(Tpre / dt), min_periods=1)

    # grab events (successful stops = true positives for training)
    print 'Fetching events (true positives)...'
    evt = fetch(dbname, 'events', *dtup[:2])['banked'].dropna()
    evt = np.around(evt / dt) * dt  # round to nearest dt
    # extend with nearby times
    truepos = (pd.DataFrame(pd.concat([evt, evt - 0.5, evt - 1.0]), 
        columns=['time']))
    truepos['outcome'] = 1

    # grab random timepoints (true negatives in training set)
    print 'Generating true negatives...'
    maxT = np.max(groupdata.index.values)
    # make some candidate random times
    Nrand = 3000
    tcands = np.random.rand(Nrand) * (maxT - Tpre) + Tpre
    tcands = np.around(tcands / dt) * dt  # round to nearest dt
    tcands = np.unique(tcands)
    # now remove all those within [Tpre, Tpost] of true positives
    rand_times = np.zeros_like(tcands)
    for cand in enumerate(tcands):
        # how close are we to already selected times?
        tt = cand[1]
        dist = truepos['time'] - tt

        # if dist in [Tpre,0] or [0, Tpost], reject
        if not (np.any(np.logical_and(dist < Tpre, dist > 0)) or 
            np.any(np.logical_and(dist > -Tpost, dist < 0))):
            rand_times[cand[0]] = tt 
    rand_times = rand_times[rand_times != 0]
    # if we have too many, we can always throw some away
    if rand_times.size > 1000:
        np.random.shuffle(rand_times)
        rand_times = rand_times[:1000]
    trueneg = pd.DataFrame(rand_times, columns=['time'])
    trueneg['outcome'] = 0

    # concatenate all training events
    allevt = pd.concat([truepos, trueneg])
    allevt = allevt.set_index('time')

    # get running average estimate of power at each timepoint of interest
    print 'Grabbing data for each event...'
    meanpwr = pd.rolling_mean(groupdata, np.ceil(Tpre / dt), min_periods=1)
    tset = pd.concat([allevt, meanpwr], axis=1, join='inner')
    tset = tset.dropna()  # can't send glmnet any row with a NaN

    # write out
    print 'Writing out...'
    outdir = '/home/jmp33/data/bartc/'
    outfile = outdir + str(dtup[0]) + '.' + str(dtup[1]) + '.lfpglmdata.csv'

    tset.to_csv(outfile)