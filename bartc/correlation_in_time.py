import numpy as np
import pandas as pd
from physutils import *
from physclasses import *
import os

# set a random seed
np.random.seed(12345)

# name of database to use
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/lfplist')

# group by (patient, dataset) entries:
groups = setlist.groupby(['patient', 'dataset'])

def get_analytic_signal(lfp, filters):
    return LFPset(lfp.demean().bandlimit(filters).decimate(5).apply(hilbert, 
        raw=True), lfp.meta.copy())

def make_correlation_series(z):
    # z is assumed to be a matrix with a channel in each column
    # return matrix is (time, chan, chan)
    zbar = np.expand_dims(np.conj(z), 2)
    return np.abs(zbar * np.expand_dims(z, 1))

def make_correlation_frame(df):
    newvals = make_correlation_series(df.values)
    # now reshape to put each pair in its own column
    newvals = np.reshape(newvals, (newvals.shape[0], -1), order='F')
    return LFPset(pd.DataFrame(newvals, index=df.index), df.meta.copy())

def get_eigens(X):
    # assume X is a vector
    Y = np.reshape(X.values, (np.sqrt(X.shape[0]), -1))
    eigs = np.linalg.eigvalsh(Y)
    return pd.Series(np.sort(eigs)[::-1])

def get_eigen_series(df):
    newdf = df.apply(get_eigens, axis=1)
    return LFPset(newdf, df.meta.copy())

# for name, grp in groups:
name = (17, 2)
lfp = fetch_all_such_LFP(dbname, *name)

# option 1: using hilbert transform and smoothing
lfpz = get_analytic_signal(lfp, ['theta'])
S = make_correlation_frame(lfpz)
Sbar = S.smooth(1).decimate(10)
eigs = get_eigen_series(Sbar) 
eigs = eigs.censor()

# option 2: doing rolling correlation, then PCA
dslfp = lfp.decimate(5)
aa = pd.rolling_corr(dslfp.dataframe, 100)