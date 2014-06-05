"""
Given a tuple of subject, dataset, line plot peri-event power in each
LFP channel
"""
from physutils import *
from physclasses import *
import rpy2.robjects as robjects
import pandas.rpy.common as com

import os
os.chdir('/home/jmp33/code/hephys/bartc')
################ all channel traces ##########################
dtup = 18, 1 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'
dt = 1./1000

# get lfp data
df = fetch_all_such(dbname, 'lfp', *dtup).set_index(['time', 'channel'])
df = df['voltage']
df = df.unstack()
    
# decimate to 40 Hz effective sampling
decfrac = (5, 5) 
dt *= np.product(decfrac)
df = dfdecimate(df, decfrac) 

# bandpass filter
df = dfbandlimit(df, ['theta'])

# instantaneous power
df = df.apply(ssig.hilbert).apply(np.absolute) ** 2

# remove censored regions
excludes = get_censor(dbname, df.index.values.astype('float64'), *dtup)
if not excludes.empty:
    excl_vec = np.any(excludes.values, axis=1)
    df[excl_vec] = np.nan

# moving average
Twin = 0.1  # size of window in s
wid = np.around(Twin / dt)
df = pd.rolling_mean(df, wid, min_periods=1, center=True)
df = df.apply(pd.Series.interpolate)

# get events
evt = fetch(dbname, 'events', *dtup)
t_evt = evt[['stop inflating', 'banked']].dropna()

# define pre and post time intervals
Tpre = -2
Tpost = 0

# grab data
chunks = df.apply(evtsplit, args=(t_evt['stop inflating'], Tpre, Tpost))

# get mean of each psth and pivot table
chanmeans = chunks.apply(pd.DataFrame.mean, axis=1)
chanmeans = chanmeans.stack().unstack(level=0)

# zscore
zscore = lambda x: (x - x.mean()) / x.std()
chanmeans = chanmeans.apply(zscore)

# prepare dataframe for passing to R
rdat = chanmeans
rdat.columns = range(rdat.shape[1])
rdat = chanmeans.reset_index()
rdat = pd.melt(rdat, id_vars='time')
rdf = com.convert_to_r_dataframe(rdat)

# load up R
R = robjects.r
R("""source('helpers.R')""")

# make plot
R("""pdf(file='~/Dropbox/hephys/media/figs/chanplot.pdf', paper='USr', width=11, height=8.5)""")

pp = R['chanmeans'](rdf)
R.plot(pp)

R('dev.off()')