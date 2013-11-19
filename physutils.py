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

def binspikes(df, dt):
    """
    Convert df, a Pandas dataframe of spike timestamps, to a binned
    histogram with bin width dt.
    """
    maxT = np.max(df['time'])
    maxbin = np.ceil(maxT / dt)
    binT = np.arange(0, maxbin * dt, dt)
    binned = np.histogram(df['time'], binT)[0]
    dframe = pd.DataFrame(binned, index=pd.Index(binT[:-1], name='time'),
        columns=['count'])
    return dframe

def smooth(df, wid):
    """
    performs smoothing by a window of width wid (in s); data in df
    reflect data at both ends to minimize edge effect
    smooths using a centered window, which is non-causal
    """
    ts = df.index[1] - df.index[0]
    x = df.values.squeeze()
    wlen = np.round(wid/ts)
    ww = np.hanning(wlen)
    # grab first wlen samples, reverse them, append to front,
    # grab last wlen samples, reverse, append to end
    xx = np.r_[x[wlen-1:0:-1], x, x[-1:-wlen:-1]]
    y = np.convolve(ww/ww.sum(),xx, mode='valid')
    y = y[(wlen/2 - 1):-wlen/2]
    return pd.DataFrame(y, index=df.index, columns=[''])

def evtsplit(df, ts, startT, endT, t0=0):
    """
    Split time series data into peri-event chunks. Data are in df.
    Times of events around which to split are in ts. 
    Code grabs startT:endT bins relative to event, so times before 
    the event have startT < 0. The first time bin is take to have 
    timestamp t0.
    """
    dt = df.index[1] - df.index[0]
    xx = df.values.squeeze()

    nevt = ts.size
    nstart = np.ceil(startT / dt)
    nend = np.ceil(endT / dt)
    binT = np.arange(nstart, nend) * dt

    evtrel = ts - t0

    elist = []
    for time in evtrel:
        offset = np.around(time / dt)
        ss = slice(nstart + offset, nend + offset)
        elist.append(pd.DataFrame(xx[ss], columns=[time]))
    alltrials = pd.concat(elist, axis=1)
    alltrials = alltrials.set_index(binT)
    alltrials.index.name = 'time'
    alltrials.columns = pd.Index(np.arange(1, nevt + 1), name='trial')
    return alltrials

def bandlimit(df, band=(0.01, 120)):
    """
    Computes bandpass-filtered version of time series in df.
    Band is either a two-element indexed sequence or a conventionally
    defined electrophysiological frequency band.
    WARNING: this is a causal filter, which incurs a phase lag!
    """
    dt = df.index[1] - df.index[0]
    band_dict = {'delta': (0.1, 4), 'theta': (4, 8), 'alpha': (8, 13), 
    'beta': (13, 30), 'gamma': (30, 100)}

    # if band isn't a two-element sequence, it should be a string
    if isinstance(band, str):
        fband = band_dict[band]

    b, a = ssig.ellip(2, 0.1, 40, [2 * dt * f for f in fband])
    bp = ssig.lfilter(b, a, df.values, axis=0)
    return pd.DataFrame(bp, index=df.index, columns=[''])

def getSpikes(*args):
    """
    Convenience function to retrieve spikes from database.
    args = (patient, dataset, channel, unit)
    Missing entries imply "return all such valid."
    """
    names = ['dataset', 'channel', 'unit']
    qstr = 'SELECT * FROM spikes WHERE patient = {}'.format(args[0])
    for tup in enumerate(args[1:]):
        qstr += ' AND {} = {}'.format(names[tup[0]], tup[1])
    qstr += ';'

    return QueryDB(qstr)

def getLFP(*args):
    """
    Convenience function to retrieve spikes from database.
    args = (patient, dataset, channel, unit)
    Missing entries imply "return all such valid."
    """
    names = ['dataset', 'channel']
    qstr = 'SELECT * FROM lfp WHERE patient = {}'.format(args[0])
    for tup in enumerate(args[1:]):
        qstr += ' AND {} = {}'.format(names[tup[0]], tup[1])
    qstr += ';'

    return QueryDB(qstr)

def getCensor(*args):
    """
    Convenience function for retrieving censoring intervals from the db
    and converting to logical arrays, one entry for each time point.
    args = patient, dataset, channel
    Assumes timestamp range equal to that of lfp.
    """
    # first, construct time axis for data by getting last two timestamps of
    # relevant lfp; use to find maximum timestamp and dt
    qstr = ("""
        SELECT time FROM lfp WHERE patient = {} AND dataset = {}
        AND channel = {} ORDER BY time DESC LIMIT 2
        """.format(*args))
    lasttwo = QueryDB(qstr).values
    dt = lasttwo[0] - lasttwo[1]
    taxis = np.arange(0, lasttwo[0] + dt, dt)

    # next, get a list of censoring times, supplement with 0 and inf 
    # on either end
    qstr = ("""
        SELECT start, stop FROM censor WHERE patient = {} 
        AND dataset = {} AND channel = {}
        """.format(*args))
    censbins = QueryDB(qstr).values.ravel()
    # supplement bins
    censbins = np.append([0.], censbins)
    censbins = np.append(censbins, np.inf)

    # finally, histogram the time axis; censored bins will have even indices
    binnum = np.digitize(taxis, censbins)
    return binnum % 2 == 0

def getEvent(event, *args):
    """
    Convenience function for retrieving task events from the db.
    event is a string of event type to retrieve
    args = patient, dataset, channel
    Assumes timestamp range equal to that of lfp.
    """

    # get events
    qstr = ("""
        SELECT {0} FROM events WHERE patient = {1} 
        AND dataset = {2} AND {0} is not NULL
        """.format(event, *args))
    return QueryDB(qstr)[event]
    


if __name__ == '__main__':

    # get all spikes for a given unit 
    df = getSpikes(18, 1, 1, 1)

    binsize = 0.050  # 50 ms bin
    binned = binspikes(df, binsize)

    evt = getEvent('banked', 18, 1)

    psth = evtsplit(binned, evt, -1, 1).mean(axis=1)

    smpsth = smooth(psth, 0.4)

    df = getLFP(18, 1, 17)
    df = df.set_index('time')