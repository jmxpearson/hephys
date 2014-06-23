import numpy as np
import pandas as pd
from physutils import *
from physclasses import *

# set a random seed
np.random.seed(12345)

# name of database to use
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/lfplist')

# group by (patient, dataset) entries:
groups = setlist.groupby(['patient', 'dataset'])

def process_LFP(lfp, filters):
    return lfp.demean().bandlimit(filters).decimate(5).apply(hilbert, 
        raw=True)

# for name, grp in groups:
name = (17, 2)
lfp = fetch_all_such_LFP(dbname, *name)
lfpz = process_LFP(lfp, ['theta'])
