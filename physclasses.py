"""
LFPset.py
LFPset is a wrapper class for a Pandas dataframe object containing LFP data.
Implements methods from physutils.py.

Methods:
decimate
bandlimit
instantaneous power
censor
split

As per Pandas convention, these should return new dataframes.
"""

import physutils
import numpy as np
from scipy.signal import hilbert
import pandas as pd

class LFPset(object):
    def __init__(self, dataframe, meta=None):
        self.dataframe = dataframe
        self.meta = meta  # dict of metadata 

    def __getattr__(self, name):
        return getattr(self.dataframe, name)

    def __str__(self):
        return self.dataframe.__str__()

    def __repr__(self):
        return 'LFP dataset object containing a\n' + self.dataframe.__repr__()

    def decimate(self, decfrac):
        newdf = physutils.dfdecimate(self.dataframe, decfrac)
        self.meta['sr'] = self.meta.get('sr', None) / np.product(decfrac)
        return LFPset(newdf, self.meta)

    def bandlimit(self, *args):
        newdf = physutils.dfbandlimit(self.dataframe, *args)
        return LFPset(newdf, self.meta)

    def demean(self):
        dmn = lambda x: (x - x.mean())
        newdf = self.dataframe.apply(dmn)
        return LFPset(newdf, self.meta)

    def zscore(self):
        zsc = lambda x: (x - x.mean()) / x.std()
        newdf = self.dataframe.apply(zsc)
        return LFPset(newdf, self.meta)

    def instpwr(self):
        newdf = self.dataframe.apply(hilbert, raw=True)
        newdf = newdf.apply(np.absolute) ** 2
        return LFPset(newdf, self.meta)

    def smooth(self, winlen):
        wid = np.around(winlen * self.meta['sr'])
        newdf = pd.rolling_mean(self.dataframe, wid, 
            min_periods=1, center=True)
        newdf = newdf.apply(pd.Series.interpolate)
        return LFPset(newdf, self.meta)

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
            return LFPset(newdf, self.meta)
        else:
            return self

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
