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
from scipy.signal import hilbert
from numpy import absolute

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
        self.meta['sr'] = self.meta.get('sr', None) / decfrac
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
        newdf = newdf.apply(absolute) ** 2
        return LFPset(newdf, self.meta)

def fetch_LFP(dbname, *tup):
    """ 
    Given a database and a tuple (tup), return an LFPset object.
    """

    lfp = physutils.fetch(dbname, 'lfp', *tup)
    lfp = lfp.set_index(['time', 'channel'])
    lfp = lfp['voltage']
    lfp = lfp.unstack()
    # the following is a kludge because the dtype is set to 'O' by the multi-index
    lfp.index = lfp.index.values.astype('float64')
    lfp.index.name = 'time'

    dt = lfp.index[1] - lfp.index[0]
    sr = 1. / dt
    meta = {'tuple': tup, 'sr': sr}

    return LFPset(lfp, meta)    

