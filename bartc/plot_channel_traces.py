"""
Given a tuple of subject, dataset, line plot peri-event power in each
LFP channel
"""
from physutils import *
from physclasses import *
import rpy2.robjects as robjects
import pandas.rpy.common as com

import os
os.chdir(os.path.expanduser('~/code/hephys/bartc'))
################ all channel traces ##########################
dtup = 18, 1 

# open data file
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# get lfp data
print "Fetching Data..."
lfp = fetch_all_such_LFP(dbname, *dtup)
    
# bandpass filter
print "Filtering..."
lfp = lfp.bandlimit(['theta'])

# decimate to 40 Hz effective sampling
print "Decimating..."
lfp = lfp.decimate(5)

# instantaneous power
print "Calculating Power..."
lfp = lfp.instpwr()

# remove censored regions
print "Censoring..."
lfp = lfp.censor()

# get events
evt = fetch(dbname, 'events', *dtup)
t_evt = evt[['stop inflating', 'banked']].dropna()

# define pre and post time intervals
Tpre = -2
Tpost = 0

# split lfp around stops
split_lfp = lfp.evtsplit(t_evt['stop inflating'], Tpre, Tpost)

# group by time and get mean for each channel for each time
chanmeans = LFPset(split_lfp.groupby(level=1).median(), meta=lfp.meta.copy()).zscore() 

print "Smoothing..."
chanmeans = chanmeans.smooth(0.2)

print "Sending to R..."
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