"""
gist to plot cross-channel LFP correlations as a function of time 
"""
from physutils import *
from physclasses import * 

# pick a dataset
dtup = 17, 2 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# get lfp data
df = fetch_all_such_LFP(dbname, *dtup)
df = df.decimate(5)
    
# step through time, take chunks, compute cross-correlation matrix