import numpy as np
import pandas as pd
import scipy.signal as ssig
from physutils import *

# define some useful numbers
np.random.seed(12345)

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/lfplist')[['patient', 'dataset']]
setlist = setlist.drop_duplicates()
dirlist = pdtbl.HDFStore(dbname).keys()

# iterate over entries:
for rec in setlist.iterrows():

    # get data indices
    dtup = rec[1].values
    # dtup = setlist.iloc[14,:].values
    print dtup

    # read in data
    print 'Reading LFP...'
    dt = 1./1000  # sampling rate in dataset
    lfp = fetch_all_such(dbname, 'lfp', *dtup, keys=dirlist)
    lfp = lfp.set_index(['time', 'channel'])
    lfp = lfp['voltage']
    lfp = lfp.unstack()
    nchan = lfp.shape[1]
    # the following is a kludge because the dtype is set to 'O' by the multi-index
    tindex = lfp.index.values.astype('float64')

    # break out by frequency band, but first de-mean
    demean = lambda x: (x - x.mean())
    lfp = lfp.apply(demean)
    print 'Filtering by frequency...'
    filters = ['delta', 'theta', 'alpha']
    allbands = dfbandlimit(lfp, filters)

    # decimate data
    print 'Decimating...'
    decfrac = (5, 5)  # reduce from 200 Hz to 40 Hz sampling rate
    for dfrac in decfrac:
        dt = dt * dfrac
        allbands = dfdecimate(allbands, dfrac)
    
    # get instantaneous power
    print 'Calculating power...'
    # could do the following with DataFrame.apply, but given how long
    # it can take, preferred option is to use the print statement to
    # show progress
    parts = []
    for aa in allbands.iteritems():
        if nchan > 1:
            print aa[0]
        parts.append(pd.DataFrame(ssig.hilbert(aa[1]), columns=[aa[0]]))
    allbands = pd.concat(parts, axis=1)
    allbands.index = tindex
    allbands.index.name = 'time'
    allbands = allbands.apply(np.absolute) ** 2

    # handle censoring
    print 'Censoring...'
    excludes = get_censor(dbname, tindex, *dtup)
    if not excludes.empty:
        excludes = excludes[excludes.columns.intersection(lfp.columns)]
        # can do something fancy later, but for now, take logical OR across all
        # channels to determine what we keep
        excl_vec = np.any(excludes.values, axis=1)
        allbands[excl_vec] = np.nan

    # standardize per channel
    print 'Standardizing regressors...'
    zscore = lambda x: (x - x.mean()) / x.std()
    allbands = allbands.apply(zscore)

    # specify peri-event times
    Tpre = 2  # time relative to event to start
    Tpost = 1.5  # time following event to exclude

    # grab events (successful stops = true positives for training)
    print 'Fetching events...'
    evt = fetch(dbname, 'events', *dtup)['banked'].dropna()
    evt = np.around(evt / dt) * dt  # round to nearest dt
    # extend with nearby times
    truepos = (pd.DataFrame(pd.concat([evt, evt - 0.5, evt - 1.0]), 
        columns=['time']))
    truepos['outcome'] = 1

    # grab random timepoints (true negatives in training set)
    print 'Generating true negatives...'
    maxT = np.max(allbands.index.values)
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
    print 'Running mean...'
    meanpwr = pd.rolling_mean(allbands, np.ceil(Tpre / dt), min_periods=1)
    tset = pd.concat([allevt, meanpwr], axis=1, join='inner')
    tset = tset.dropna()  # can't send glmnet any row with a NaN

    # write out
    print 'Writing out...'
    outdir = '/home/jmp33/data/bartc/'
    outfile = outdir + str(dtup[0]) + '.' + str(dtup[1]) + '.lfpglmdata.csv'

    tset.to_csv(outfile)