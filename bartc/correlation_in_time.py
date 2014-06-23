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

def get_analytic_signal(lfp, filters):
    return lfp.demean().bandlimit(filters).decimate(5).apply(hilbert, 
        raw=True)

def make_correlation_matrix(z):
    # z is assumed to be a matrix with a channel in each column
    # return matrix is (time, chan, chan)
    zbar = np.expand_dims(np.conj(z), 2)
    return np.abs(zbar * np.expand_dims(z, 1))

# for name, grp in groups:
name = (17, 2)
lfp = fetch_all_such_LFP(dbname, *name)
lfpz = get_analytic_signal(lfp, ['theta']).values
S = make_correlation_matrix(lfpz)
