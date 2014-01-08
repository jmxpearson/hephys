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

class LFPset(object):
    def __init__(self, dataframe, sr=1000):
        self.dataframe = dataframe
        self.sr = sr  # sampling rate in Hz

    def __getattr__(self, name):
        return getattr(self.dataframe, name)

    def __str__(self):
        return self.dataframe.__str__()

    def __repr__(self):
        return 'LFP dataset object containing a\n' + self.dataframe.__repr__()

    def decimate(self, decfrac):
        newdf = physutils.dfdecimate(self.dataframe, decfrac)
        return LFPset(newdf, self.sr / decfrac)

    def bandlimit(self, *args):
        newdf = physutils.bandlimit(self.dataframe, *args)
        return LFPset(newdf, self.sr)
