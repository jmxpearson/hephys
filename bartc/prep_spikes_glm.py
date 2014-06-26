"""
For each patient, dataset, unit, prepare a dataframe to be used to fit a 
generalized linear model.
"""
import numpy as np
import os
import pandas as pd
from physutils import *

def make_regressor_is_in_trial(taxis, events):
    # takes a dataframe of events and a time axis and returns a binary series
    is_in_trial = pd.Series(0, index=taxis, name='is_in_trial')
    starts = evt['trial_start']
    stops = evt['trial_over']
    pairs = zip(starts, stops)
    for p in pairs:
        is_in_trial[slice(*p)] = 1
    return is_in_trial


# set a random seed
np.random.seed(12345)

# name of database to use
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/spklist')

# for idx, row in setlist.iterrows():
dtup = (17, 2, 1, 1)

spks = load_spikes(dbname, dtup)

evt = fetch(dbname, 'events', *dtup[0:2])

regressor_list = [make_regressor_is_in_trial]

regressor_frame = pd.concat([f(spks.index, evt) for f in regressor_list], 
    axis=1)