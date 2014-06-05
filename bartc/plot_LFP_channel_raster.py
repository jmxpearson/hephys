"""
Contains plotting functions, some of which read data from the 
database and call R for plotting.
"""
from physutils import *
from physclasses import *
import rpy2.robjects as robjects
import pandas.rpy.common as com

import os
os.chdir('/home/jmp33/code/hephys/bartc')

################ single channel raster ##########################
dtup = 18, 1, 22 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# get lfp data
lfpset = fetch_LFP(dbname, *dtup)

# bandpass filter
lfpset = lfpset.bandlimit(['theta'])

# decimate to 40 Hz effective sampling
lfpset = lfpset.decimate(5)

# instantaneous power
lfpset = lfpset.instpwr()

# remove censored regions
lfpset = lfpset.censor()

# rolling smooth
lfpset = lfpset.smooth(0.1)

# convert to dB
lfpset = lfpset.apply(lambda x: 10 * np.log10(x))

# get events
evt = fetch(dbname, 'events', *dtup[:-1])
max_inf = evt['inflate_time'].max()
t_evt = evt[['start inflating', 'inflate_time','banked']].dropna()

# define pre and post time intervals
Tpre = -1
Tpost = np.ceil(max_inf) + 1

# grab data
chunks = evtsplit(lfpset, t_evt['start inflating'], Tpre, Tpost)

# sort trials by inflate time
sortord = np.argsort(t_evt['inflate_time']) + 1
chunks = chunks[sortord].reset_index()
chunks = pd.melt(chunks, id_vars=['time'])

# make a dataframe of inflate times for plot annotation
inflate_times = t_evt['inflate_time'].copy()
inflate_times.sort()
inflate_times = inflate_times.reset_index(drop=True)
inflate_times = pd.DataFrame(inflate_times).reset_index()

# load up R
R = robjects.r
R('source("helpers.R")')
rdf = com.convert_to_r_dataframe(chunks)
inft = com.convert_to_r_dataframe(inflate_times)

# make plot
R("""pdf(file='~/Dropbox/hephys/media/figs/powraster.pdf', paper='USr', width=11, height=8.5)""")

pp = R['rasterize'](rdf, inft)
R.plot(pp)

R('dev.off()')

