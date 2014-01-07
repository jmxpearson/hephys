"""
gist to plot cross-channel LFP correlations as a function of time 
"""
from physutils import * 

# pick a dataset
dtup = 17, 2 

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