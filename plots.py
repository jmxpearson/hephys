"""
Contains plotting functions, some of which read data from the 
database and call R for plotting.
"""
from physutils import *
import rpy2.robjects as robjects
import pandas.rpy.common as com

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'
dt = 1./200

# let's try to make a raster dataframe
dtup = 18, 1, 22 

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