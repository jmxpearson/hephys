"""
Given a tuple of subject, dataset, line plot peri-event power in each
LFP channel
"""
import rpy2.robjects as robjects
import pandas.rpy.common as com
import pandas as pd
import physutils
from hephys import dbio
import os

os.chdir(os.path.expanduser('~/code/hephys/bartc'))
################ all channel traces ##########################
dtup = 20, 1 

# open data file
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# get lfp data
print "Fetching Data..."
lfp = dbio.fetch_all_such_LFP(dbname, *dtup)
    
# bandpass filter
print "Filtering..."
lfp = lfp.bandlimit(['theta'])

# decimate to 100 Hz effective sampling
print "Decimating..."
lfp = lfp.decimate(5)

# instantaneous power
print "Calculating Power..."
lfp = lfp.instpwr()

# remove censored regions
print "Censoring..."
lfp = lfp.censor()

# get events
evt = dbio.fetch(dbname, 'events', *dtup)
t_evt = evt[['stop inflating', 'banked']].dropna()
stops = t_evt['stop inflating']
pops = evt['popped'].dropna()

# define pre and post time intervals
Tpre = -2
Tpost = 0

# split lfp around stops
stop_split = lfp.evtsplit(stops, Tpre, Tpost)
pop_split = lfp.evtsplit(pops, Tpre, Tpost)

# group by time and get mean for each channel for each time
stop_means = physutils.LFPset(stop_split.groupby(level=1).median(), meta=lfp.meta.copy()).zscore() 
pop_means = physutils.LFPset(pop_split.groupby(level=1).median(), meta=lfp.meta.copy()).zscore() 

print "Smoothing..."
stop_means= stop_means.smooth(0.2)
pop_means= pop_means.smooth(0.2)

print "Sending to R..."
def make_filename(df, name):
    base = '~/Dropbox/hephys/media/figs/' 
    colnames = df.columns
    bands = set(map(lambda x: x.split('.')[0], colnames))
    pieces = list(df.meta['tuple']) + ['_'.join(bands), name, 'chanplot', 'pdf']
    return "\'" + base + '.'.join(map(str, pieces)) + "\'"

def print_from_R(df, name):
    # prepare dataframe for passing to R
    rdat = df.copy()
    rdat.columns = range(rdat.shape[1])
    rdat = rdat.reset_index()
    rdat = pd.melt(rdat, id_vars='time')
    rdf = com.convert_to_r_dataframe(rdat)

    fname = make_filename(df, name)

    # load up R
    R = robjects.r
    R("""source('helpers.R')""")

    # make plot
    R("pdf(file=" + fname + ", paper='USr', width=11, height=8.5)")

    pp = R['chanmeans'](rdf)
    R.plot(pp)

    R('dev.off()')

frames = [stop_means, pop_means]
names = ['stop', 'pop']
for (frame, name) in zip(frames, names):
    print_from_R(frame, name)