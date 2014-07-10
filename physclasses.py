"""
LFPset.py
LFPset is a wrapper class for a Pandas dataframe object containing LFP data.
Implements methods from physutils.py.

Methods:
decimate
bandlimit
instantaneous power
censor
evtsplit

As per Pandas convention, these should return new dataframes.
"""

import physutils
import numpy as np
from scipy.signal import hilbert
import pandas as pd

class LFPset(object):
    def __init__(self, data, meta=None):
        # wrap passed data in constructor in case it's a series
        self.dataframe = pd.DataFrame(data)
        self.meta = meta  # dict of metadata 

    def __getattr__(self, name):
        return getattr(self.dataframe, name)

    def __str__(self):
        return self.dataframe.__str__()

    def __repr__(self):
        return 'LFP dataset object containing a\n' + self.dataframe.__repr__()

    def decimate(self, decfrac):
        newdf = physutils.dfdecimate(self.dataframe, decfrac)
        newmeta = self.meta.copy()
        newmeta['sr'] = newmeta.get('sr', None) / np.product(decfrac)
        return LFPset(newdf, newmeta)

    def bandlimit(self, *args):
        newdf = physutils.dfbandlimit(self.dataframe, *args)
        newmeta = self.meta.copy()
        return LFPset(newdf, newmeta)

    def demean(self):
        dmn = lambda x: (x - x.mean())
        newdf = self.dataframe.apply(dmn)
        newmeta = self.meta.copy()
        return LFPset(newdf, newmeta)

    def interpolate(self):
        newdf = self.dataframe.interpolate()
        newmeta = self.meta.copy()
        return LFPset(newdf, newmeta)

    def zscore(self):
        zsc = lambda x: (x - x.mean()) / x.std()
        newdf = self.dataframe.apply(zsc)
        newmeta = self.meta.copy()
        return LFPset(newdf, newmeta)

    def instpwr(self):
        # by padding up to next highest power of 2, we get a huge 
        # performance boost; truncate afterward
        Nstart = self.dataframe.shape[0]
        Nfft = 2 ** np.ceil(np.log2(Nstart))
        hilbert_pad = lambda x: hilbert(x, N=Nfft)[:Nstart]
        newdf = self.dataframe.apply(hilbert_pad, raw=True)
        newdf = newdf.apply(np.absolute) ** 2
        newmeta = self.meta.copy()
        return LFPset(newdf, newmeta)

    def smooth(self, winlen):
        wid = np.around(winlen * self.meta['sr'])
        newdf = pd.rolling_mean(self.dataframe, wid, 
            min_periods=1)
        newdf = newdf.apply(pd.Series.interpolate)
        newmeta = self.meta.copy()
        return LFPset(newdf, newmeta)

    def censor(self):
        excludes = physutils.get_censor(
            self.meta['dbname'], self.dataframe.index, *self.meta['tuple'])
        if not excludes.empty:
            excludes = excludes[excludes.columns.intersection(
                self.dataframe.columns)]
            # can do something fancy later, but for now, take logical OR across all
            # channels to determine what we keep
            excl_vec = np.any(excludes.values, axis=1)
            newdf = self.dataframe
            newdf[excl_vec] = np.nan
            newmeta = self.meta.copy()
            return LFPset(newdf, newmeta)
        else:
            return self

    def evtsplit(self, times, Tpre, Tpost, t0=0):
        # note: Tpre < 0 for times before the events in times
        split_to_series = lambda x: physutils.evtsplit(x, times, Tpre, Tpost, 
            t0).unstack()
        return self.dataframe.apply(split_to_series)

    def avg_time_frequency(self, channel, times, Tpre, Tpost, method='wav', doplot=True, **kwargs):
        """
        Do a time frequency decomposition, averaging across chunks split at 
        times. **kwargs are passed on as parameters to the method call:
        method 'wav': w parameter
        method 'spec': window length (in s) and fraction overlap
        """
        series = self.dataframe[channel]
        if method == 'wav':
            callback = physutils.continuous_wavelet
        else:
            callback = physutils.spectrogram

        tf = physutils.avg_time_frequency(series, callback, times, Tpre, 
            Tpost, **kwargs)

        if doplot:
            fig = physutils.plot_time_frequency(tf)
        else:
            fig = None

        return tf, fig



def fetch_LFP(dbname, *tup):
    """ 
    Given a database and a tuple (tup), return an LFPset object.
    """

    lfp = physutils.fetch(dbname, 'lfp', *tup)
    lfp = lfp.set_index(['time', 'channel'])
    lfp = lfp['voltage']
    lfp = lfp.unstack()
    
    dt = lfp.index[1] - lfp.index[0]
    sr = 1. / dt
    meta = {'dbname': dbname, 'tuple': tup, 'sr': sr}

    return LFPset(lfp, meta)    

def fetch_all_such_LFP(dbname, *tup, **kwargs):
    """ 
    Given a database and a tuple (tup), return an LFPset object.
    """

    lfp = physutils.fetch_all_such(dbname, 'lfp', *tup, **kwargs)
    lfp = lfp.set_index(['time', 'channel'])
    lfp = lfp['voltage']
    lfp = lfp.unstack()
    
    dt = lfp.index[1] - lfp.index[0]
    sr = 1. / dt
    meta = {'dbname': dbname, 'tuple': tup, 'sr': sr}

    return LFPset(lfp, meta)
