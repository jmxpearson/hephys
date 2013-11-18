import numpy as np
import MySQLdb
import pandas as pd
import pandas.io.sql as pdsql
import scipy.signal as ssig
import warnings

def decimate(x, decfrac, axis=-1):
    """
    just like in the scipy source code, except I use filtfilt
    and opt for hardcode defaults
    q is the fraction to decimate by; length of returned data is len(x)/q
    """
    if decfrac > 15:
        wrnstr = """You are attempting to decimate by a factor > 15. You are risking numerical instability. Consider performing multiple successive decimations instead."""
        warnings.warn(wrnstr)

    n = 8
    b, a = ssig.filter_design.cheby1(n, 0.05, 0.8/decfrac)

    y = ssig.filtfilt(b, a, x, axis=axis)

    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, decfrac)
    return y[sl]

def QueryDB(qstr):
    """
    Query database using SQL query in qstr.
    """
    db = MySQLdb.connect(host='localhost',
        user='root', passwd='', db='bartc')
    dat = pdsql.read_frame(qstr, db)
    return dat

def binarize(df, dt):
    """
    Convert df, a Pandas dataframe of spike timestamps, to a binned
    histogram with bin width dt.
    """
    maxT = np.max(df['time'])
    maxbin = np.ceil(maxT / dt)
    binT = np.arange(0, maxbin * dt, dt)
    binned = np.histogram(df['time'], binT)[0]
    return pd.DataFrame(binned, index=binT[:-1])




# get all spikes for a given unit 
qstr = """SELECT *
FROM spikes WHERE patient = 18 AND dataset = 1 AND channel = 1 AND unit = 1"""

binsize = 0.050  # 50 ms bin