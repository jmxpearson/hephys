import numpy as np
import MySQLdb
import pandas as pd
import pandas.io.sql as pdsql
import scipy.signal as ssig
from physutils import *

# define some useful numbers
dt = 1./200  # sampling rate in dataset

# first, get a list of lfp channels
qstr = """SELECT DISTINCT patient, dataset, channel FROM lfp;"""
chanlist = QueryDB(qstr)



# get data indices
dtup = chanlist.iloc[49,:].values

# read in data
lfp = getLFP(*dtup)
lfp = lfp.set_index('time')
lfp = lfp['voltage']

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
allbands.columns = filters

# handle censoring
excludes = getCensor(*dtup)
excludes = excludes[::decfrac]
allbands[excludes] = np.nan

# grab events (successful stops: true positives for training)
evt = getEvent('banked', *dtup[:-1])
# extend with nearby times
evt = pd.DataFrame(pd.concat([evt, evt - 0.5, evt - 1.0]), columns=['time'])
evt['outcome'] = 1

# specify peri-event times
Tpre = -2  # time relative to event to start
Tpost = 1.5  # time following event to exclude

# grab random timepoints (true negatives in training set)
Nrand = 1000
maxT = np.max(allbands.index.values)
rand_times = np.array([])
while len(rand_times) < Nrand:
    # generate a candidate time point with at least abs(Tpre) before to grab
    tt = np.random.rand() * (maxT + Tpre) - Tpre

    # how close are we to already selected times
    dist = evt - tt
    #rdist = np.abs(rand_times - tt)

    # if dist in [Tpre,0] or [0, Tpost], reject
    if (np.any(np.logical_and(dist < -Tpre, dist > 0)) or 
        np.any(np.logical_and(dist > -Tpost, dist < 0))): 
        pass
    else:
        rand_times = np.append(rand_times, tt)
trueneg = pd.DataFrame(rand_times, columns=['time'])
trueneg['outcome'] = 0

# concatenate all training events
evt = pd.concat([evt, trueneg])