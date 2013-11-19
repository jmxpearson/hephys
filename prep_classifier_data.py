import numpy as np
import pandas as pd
import scipy.signal as ssig
from physutils import *

# define some useful numbers
dt = 1./200  # sampling rate in dataset

# first, get a list of lfp channels
qstr = """SELECT DISTINCT patient, dataset FROM events;"""
setlist = QueryDB(qstr)



# get data indices
dtup = setlist.iloc[14,:].values
print dtup

# read in data
lfp = getLFP(*dtup)
lfp = lfp.set_index(['time', 'channel'])
lfp = lfp['voltage']
lfp = lfp.unstack()

# break out by frequency band
filters = ['delta', 'theta', 'alpha']
bands = [bandlimit(lfp, f) for f in filters]
allbands = pd.concat(bands, axis=1)

# decimate data
decfrac = 5  # reduce from 200 Hz to 40 Hz sampling rate
dt = dt * decfrac
tindex = np.arange(0, allbands.index[-1], dt)
parts = [pd.DataFrame(decimate(aa[1], decfrac)) for aa in allbands.iteritems()]
allbands = pd.concat(parts, axis=1)
allbands.index = tindex
allbands.index.name = 'time'

# get instantaneous power
allbands = allbands.apply(ssig.hilbert)
allbands = allbands.apply(np.abs)

# handle censoring
excludes = getCensor(*dtup)
excludes = excludes[::decfrac]
allbands[excludes] = np.nan

# grab events (successful stops = true positives for training)
evt = getEvent('banked', *dtup[:-1])
# extend with nearby times
truepos = (pd.DataFrame(pd.concat([evt, evt - 0.5, evt - 1.0]), 
    columns=['time']))
truepos['outcome'] = 1

# specify peri-event times
Tpre = 2  # time relative to event to start
Tpost = 1.5  # time following event to exclude

# grab random timepoints (true negatives in training set)
Nrand = 1000
maxT = np.max(allbands.index.values)
rand_times = np.array([])
while len(rand_times) < Nrand:
    # generate a candidate time point with at least Tpre before to grab
    tt = np.random.rand() * (maxT - Tpre) + Tpre

    # how close are we to already selected times
    dist = evt - tt
    #rdist = np.abs(rand_times - tt)

    # if dist in [Tpre,0] or [0, Tpost], reject
    if (np.any(np.logical_and(dist < Tpre, dist > 0)) or 
        np.any(np.logical_and(dist > -Tpost, dist < 0))): 
        pass
    else:
        rand_times = np.append(rand_times, tt)
trueneg = pd.DataFrame(rand_times, columns=['time'])
trueneg['outcome'] = 0

# concatenate all training events
allevt = pd.concat([truepos, trueneg])
allevt['time'] = np.around(allevt['time'] / dt) * dt
allevt = allevt.set_index('time')

# get running average estimate of power at each timepoint of interest
meanpwr = pd.rolling_mean(allbands, np.ceil(Tpre / dt), min_periods=1)
tset = pd.concat([meanpwr, allevt], axis=1, join='inner')

# write out
outdir = '/home/jmp33/data/bartc/'
outfile = outdir + str(dtup[0]) + '.' + str(dtup[1]) + '.lfpglmdata.csv'

tset.to_csv(outfile)