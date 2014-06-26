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

def corr_then_med(lfp, events, Tpre, Tpost):
    # get correlation matrix between time series, extract eigenvalues, 
    # median across trials
    # Tpre and Tpost are times relative to events to grab

    print "Decimating..."
    lfpd = lfp.decimate(5)  # decimate to 40 Hz

    print "Extracting eigenvalues..."
    eigseries = PCA_series(lfpd, win=0.25, frac_overlap=0.95).censor()

    # norm to max eigenvalue
    total_var_normed = LFPset(eigseries.div(eigseries.dataframe.iloc[:, -1], axis=0), meta=eigseries.meta.copy())

    # normalize to max eigenvalue
    max_eigen_normed = LFPset(eigseries.div(eigseries.sum(axis=1), axis=0), meta=eigseries.meta.copy())

    # break up data around events
    frames = (eigseries, total_var_normed, max_eigen_normed)
    print "Splitting Data..."
    splits = [df.evtsplit(events, Tpre, Tpost) for df in frames]

    # group by time and get mean for each channel for each time
    medians = [df.groupby(level=1).median() for df in splits]
    return tuple(medians) 

def med_then_corr(lfp, events, Tpre, Tpost):
    # median channels across trials, then get correlation matrix across time
    # and extract eigenvalue series

    print "Decimating..."
    lfpd = lfp.decimate(5).censor()  # decimate to 40 Hz

    # break up data around events
    print "Splitting Data..."
    split_lfp = lfpd.evtsplit(events, Tpre, Tpost)

    chanmeds = LFPset(split_lfp.groupby(level=1).median(), 
        meta=lfpd.meta.copy())

    print "Extracting eigenvalues..."
    eigseries = PCA_series(chanmeds, win=0.25, frac_overlap=0.95)
    
    # normlize as percent total variance
    total_var_normed = LFPset(eigseries.div(eigseries.sum(axis=1), axis=0), meta=eigseries.meta.copy())

    # normalize to max eigenvalue
    max_eigen_normed = LFPset(eigseries.div(eigseries.dataframe.iloc[:, -1], axis=0), meta=eigseries.meta.copy())

    return (eigseries, total_var_normed, max_eigen_normed)

# for name, grp in groups:
name = (18, 1)
print name

print "Fetching LFP data..."
lfp = fetch_all_such_LFP(dbname, *name)

# fetch events
evt = fetch(dbname, 'events', *name)
stops = evt[['stop inflating']].dropna().values

# prepare to do some plotting
Tpre = -2
Tpost = 0

median_of_corr = corr_then_med(lfp, stops, Tpre, Tpost)
corr_of_median = med_then_corr(lfp, stops, Tpre, Tpost)
