import numpy as np
import pandas as pd
import os
from physutils import *
from physclasses import *

# set a random seed
np.random.seed(12345)

# name of database to use
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/lfplist')

# group by (patient, dataset) entries:
groups = setlist.groupby(['patient', 'dataset'])

def get_eigens(X):
    # assume X is a vector
    numeigs = np.sqrt(X.shape[0])
    Y = np.reshape(X.values, (numeigs, -1))
    if np.any(np.isnan(Y)):
        eigs = np.empty(numeigs)
        eigs[:] = np.nan
    else:
        eigs = np.linalg.eigvalsh(Y)
        eigs = np.sort(eigs[::-1])
    return pd.Series(eigs)

def get_eigen_series(df):
    newdf = df.apply(get_eigens, axis=1)
    return LFPset(newdf, df.meta.copy())

def make_PCA_series(lfp, window, skip):
    # window is the length of the window to slide (in bins)
    # skip is the spacing between samples in the returned frame
    corrseries = pd.rolling_corr(lfp.dataframe, window=window)
    corrseries = corrseries.iloc[::skip, :, :]
    corrseries = corrseries.to_frame(filter_observations=False).transpose()
    lfp_pairs = LFPset(corrseries, meta=lfp.meta.copy())
    eigs = get_eigen_series(lfp_pairs)
    return eigs

def PCA_series(lfp, win, frac_overlap=1):
    # win is the length of the running window (in s)
    # frac_overlap is the fraction of overlap between successive 
    # samples returned
    dt = 1. / lfp.meta['sr']
    winbins = int(np.around(win / dt))
    skip = max(1, int((1 - frac_overlap) * winbins))
    return make_PCA_series(lfp, winbins, skip)

# for name, grp in groups:
name = (17, 2)
lfp = fetch_all_such_LFP(dbname, *name)
lfp = lfp.decimate(5)  # decimate to 40 Hz
eigseries = PCA_series(lfp, win=2.0, frac_overlap=0.75).censor()
var_explained = eigseries.div(eigseries.sum(axis=1), axis=0)

# fetch events
evt = fetch(dbname, 'events', *name)
t_evt = evt[['stop inflating', 'banked']].dropna()

# prepare to do some plotting
Tpre = -2
Tpost = 0

# break up data around events
split_lfp = eigseries.evtsplit(t_evt['stop inflating'], Tpre, Tpost)

# # get mean of each psth and pivot table
# chanmeans = chunks.apply(pd.DataFrame.mean, axis=1)
# chanmeans = chanmeans.stack().unstack(level=0)

# # zscore
# zscore = lambda x: (x - x.mean()) / x.std()
# chanmeans = chanmeans.apply(zscore)
