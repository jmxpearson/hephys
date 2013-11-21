"""
Contains plotting functions, some of which read data from the 
database and call R for plotting.
"""
from physutils import *
import rpy2.robjects as robjects
import pandas.rpy.common as com

################ single channel raster ##########################
dtup = 18, 1, 22 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'
dt = 1./200

# get lfp data
df = fetch(dbname, 'lfp', *dtup).set_index('time')['voltage']

# censor
excludes = get_censor(dbname, df.index.values.astype('float64'), *dtup[:2])
if not excludes.empty:
    excl_vec = np.any(excludes.values, axis=1)
    
# decimate to 40 Hz effective sampling
decfrac = 5
dt *= decfrac
df = dfdecimate(df, decfrac) 

# smooth censoring region, too
excl_vec = dfdecimate(excl_vec.astype('int'), decfrac) > 0.8

# bandpass filter
df = dfbandlimit(df, ['theta'])

# instantaneous power
df = df.apply(ssig.hilbert).apply(np.absolute) ** 2

# remove censored regions
df[excl_vec] = np.nan

# rolling smooth
Twin = 0.1  # size of window in s
wid = np.around(Twin / dt)
df = pd.rolling_mean(df, wid, min_periods=1, center=True)
df = df.apply(pd.Series.interpolate)

# convert to dB
df = df.apply(lambda x: 10 * np.log10(x))

# get events
evt = fetch(dbname, 'events', *dtup[:-1])
max_inf = evt['inflate_time'].max()
t_evt = evt[['start inflating', 'inflate_time','banked']].dropna()

# define pre and post time intervals
Tpre = -1
Tpost = np.ceil(max_inf) + 1

# grab data
chunks = evtsplit(df, t_evt['start inflating'], Tpre, Tpost)

# sort trials by inflate time
sortord = np.argsort(t_evt['inflate_time']) + 1
chunks = chunks[sortord].reset_index()
chunks = pd.melt(chunks, id_vars=['time'])

# make a dataframe of inflate times for plot annotation
inflate_times = t_evt['inflate_time'].copy()
inflate_times.sort()
inflate_times = inflate_times.reset_index(drop=True)
inflate_times = pd.DataFrame(inflate_times).reset_index()
# chunks.to_csv('/home/jmp33/data/bartc/testchunks.csv')
# t_evt['inflate_time'].to_csv('/home/jmp33/data/bartc/inftimes.csv')

# load up R
R = robjects.r
rdf = com.convert_to_r_dataframe(chunks)
inft = com.convert_to_r_dataframe(inflate_times)

R('source("helpers.R")')
pp = R['rasterize'](rdf, inft)
R.X11()
R.plot(pp)

################ all channel traces ##########################
dtup = 18, 1 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'
dt = 1./200

# get lfp data
df = fetch_all_such(dbname, 'lfp', *dtup).set_index(['time', 'channel'])
df = df['voltage']
df.unstack()
    
# decimate to 40 Hz effective sampling
decfrac = 5
dt *= decfrac
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