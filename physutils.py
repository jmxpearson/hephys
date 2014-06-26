import h5py
import numpy as np
import pandas as pd
import pandas.io.pytables as pdtbl
import scipy.signal as ssig
import warnings

def make_path(*tup):
    abbr = ['p', 'd', 'c', 'u'][:len(tup)]
    nstrs = map(str, tup)
    pieces = [a + b for a,b in zip(abbr, nstrs)]
    return '/'.join(pieces)

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
    b, a = ssig.filter_design.cheby1(n, 0.05, 0.8 / decfrac)

    y = ssig.filtfilt(b, a, x, axis=axis)

    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, decfrac)
    return y[sl]

def dfdecimate(df, decfrac):
    """
    Decimate a dataframe, handling indices and columns appropriately.
    decfrac can be an iterable of successive decimations
    """
    newdf = pd.DataFrame(df)  # upcast from Series, if needed

    # if we passed an int, make it a tuple
    if isinstance(decfrac, int):
        decfrac = (decfrac,)

    for frac in decfrac:
        tindex = newdf.index[::frac]
        parts = [pd.DataFrame(decimate(aa[1], frac), columns=[aa[0]]) 
        for aa in newdf.iteritems()]
        newdf = pd.concat(parts, axis=1)
        newdf.index = tindex
        newdf.index.name = df.index.name
    return newdf

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
    the event have startT < 0. The first time bin is assumed to have 
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
    alltrials.columns = pd.Index(np.arange(nevt), name='trial')
    return alltrials

def bandlimit(df, band=(0.01, 120)):
    """
    Computes bandpass-filtered version of time series in df.
    Band is either a two-element indexed sequence or a conventionally
    defined electrophysiological frequency band.
    """
    dt = df.index[1] - df.index[0]
    band_dict = {'delta': (0.1, 4), 'theta': (4, 8), 'alpha': (8, 13), 
    'beta': (13, 30), 'gamma': (30, 100)}

    # if band isn't a two-element sequence, it should be a string
    if isinstance(band, str):
        fband = band_dict[band]
    else:
        fband = band

    b, a = ssig.iirfilter(2, [2 * dt * f for f in fband], rp=0.1, rs=40,
        ftype='ellip')
    return df.apply(lambda x: ssig.filtfilt(b, a, x), raw=True)

def dfbandlimit(df, filters=None):
    """
    Convenience function for bandlimiting data frames. Handles
    indices and columns appropriately.
    """
    if filters is None:
        return df

    df = pd.DataFrame(df)
    nchan = df.shape[1]
    bands = [bandlimit(df, f) for f in filters]
    allbands = pd.concat(bands, axis=1)
    
    # attend to labeling
    fstr = map(str, filters)
    bandpairs = zip(np.repeat(fstr, nchan), allbands.columns)
    bandnames = [b[0] + '.' + str(b[1]) for b in bandpairs]
    allbands.columns = bandnames

    return allbands

def fetch(dbname, node, *args):
    """
    Given a node ('lfp', 'spikes', 'events', 'censor'), and
    a tuple of (patient, dataset, channel, unit), retrieves data.
    """
    target = node + '/' + make_path(*args)
    return pd.read_hdf(dbname, target)

def fetch_metadata(dbname, node, *args):
    # have to get metadata this way because we didn't store it via pandas
    # (which currently has no support for dataframe metadata)
    target = node + '/' + make_path(*args)
    fobj = h5py.File(dbname, 'r')
    attdict = dict(fobj[target].attrs)
    fobj.close()
    return attdict

def fetch_all_such(dbname, node, *args, **kwargs):
    """
    Given an incomplete specification in args, get all datasets consistent
    with it. Return a single dataframe.
    If database keys are precomputed to save time, they can be specified.
    """
    if 'keys' in kwargs:
        keys = kwargs['keys']
    else:
        keys = pdtbl.HDFStore(dbname).keys()

    glob = node + '/' + make_path(*args)
    # do really simple regex matching
    matches = [k for k in keys if glob in k]

    if matches:
        parts = []
        for m in matches:
            parts.append(fetch(dbname, m))
        excl = pd.concat(parts)
    else:
        excl = pd.DataFrame([])

    return excl

def get_censor(dbname, taxis, *args):
    """
    Convenience function for retrieving censoring intervals from the db
    and converting to logical arrays, one entry for each time point.
    args = patient, dataset, channel
    Assumes timestamp range equal to that of lfp.
    """
    
    # make sure tindex is of type float64
    ##### should be able to remove when pandas is upgraded to 0.13, which 
    ##### allows for float64 indices
    taxis = pd.Series(taxis.values.astype('float64'))

    # get censoring intervals, group by channel
    censors = fetch_all_such(dbname, 'censor', *args)
    if not censors.empty:
        censors = censors.groupby('channel')
    else:
        return censors
    
    # arrange start and stop times into linear sequence
    # (the extra pair of braces around the return value is to prevent
    # pandas from converting the time array to a series when apply gets only
    # a single return value (i.e., when censors has only a single group))
    if censors.ngroups > 1:
        flatfun = lambda x: x[['start', 'stop']].values.ravel()
    else:
        flatfun = lambda x: [x[['start', 'stop']].values.ravel()]
    censbins = censors.apply(flatfun)
    
    # append 0 and inf to bins
    censbins = censbins.apply(lambda x: np.append([0], x))
    censbins = censbins.apply(lambda x: np.append(x, np.inf))
    # bin times in taxis; censored bins will have even indices
    binnum = censbins.apply(lambda x: np.digitize(taxis, x))
    binnum = binnum.apply(lambda x: x % 2 == 0)

    excludes = pd.concat([pd.Series(b) for b in binnum], axis=1)
    excludes.columns = binnum.index  # channel names
    excludes.index = taxis
    
    return excludes
   
def censor_spikes(df, dbname, dtup):
    excludes = get_censor(dbname, df.index, *dtup)
    if not excludes.empty:
        excludes = excludes[excludes.columns.intersection(
            df.columns)]
        excl_vec = np.any(excludes.values, axis=1)
        newdf = df.copy()
        newdf[excl_vec] = np.nan
        return newdf
    else:
        return df

def load_spikes(dbname, dtup):
    spks = fetch(dbname, 'spikes', *dtup)

    spkbin = binspikes(spks, 0.05)  # use 50 ms bins

    return censor_spikes(spkbin, dbname, dtup)

if __name__ == '__main__':

    # get all spikes for a given unit 
    # df = getSpikes(18, 1, 1, 1)
    dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'
    df = fetch(dbname, 'spikes', 18, 1, 1, 1)

    binsize = 0.050  # 50 ms bin
    binned = binspikes(df, binsize)

    evt = fetch(dbname, 'events', 18, 1)['banked'].dropna()

    psth = evtsplit(binned, evt, -1, 1).mean(axis=1)

    smpsth = smooth(psth, 0.4)

    df = fetch_all_such(dbname, 'spikes', 17, 2)

    df = fetch(dbname, 'lfp', 18, 1, 17)
    df = df.set_index('time')