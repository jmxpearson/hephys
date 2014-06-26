"""
For each patient, dataset, unit, prepare a dataframe to be used to fit a 
generalized linear model.
"""
import numpy as np
import os
import pandas as pd
from physutils import *

# set a random seed
np.random.seed(12345)

# name of database to use
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/spklist')

# for idx, row in setlist.iterrows()
dtup = (17, 2, 1, 1)

spks = fetch(dbname, 'spikes', *dtup)