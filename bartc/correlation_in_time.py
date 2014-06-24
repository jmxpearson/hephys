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
    corrseries = pd.rolling_corr(lfp.dataframe, window=window, 
        min_periods=1)
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


