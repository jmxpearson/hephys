"""
LFPset.py
LFPset is a wrapper class for a Pandas dataframe object containing LFP data.

Methods:
decimate
bandlimit
instantaneous power
censor
split

As per Pandas convention, these should return new dataframes.
"""

class LFPset:
    def __init__(self, dataframe):
        self.dataframe = dataframe

    def __getattr__(self, name):
        return getattr(self.dataframe, name)
