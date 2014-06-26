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
    reg = pd.Series(0, index=taxis, name='is_in_trial')
    starts = events['trial_start']
    stops = events['trial_over']
    pairs = zip(starts, stops)
    for p in pairs:
        reg[slice(*p)] = 1
    return reg

def make_regressor_is_inflating(taxis, events):
    # takes a dataframe of events and a time axis and returns a binary series
    reg = pd.Series(0, index=taxis, name='is_inflating')
    starts = events['start inflating']
    stops = events[['stop inflating', 'popped']].sum(axis=1)
    pairs = zip(starts, stops)
    for p in pairs:
        reg[slice(*p)] = 1
    return reg

def make_regressor_elapsed_time(taxis, events):
    # takes a dataframe of events and a time axis and returns a binary series
    reg = pd.Series(0, index=taxis, name='elapsed_time')
    starts = events['start inflating']
    stops = events[['stop inflating', 'popped']].sum(axis=1)
    pairs = zip(starts, stops)
    for p in pairs:
        slc = slice(*p)
        time = reg[slc].index.values
        reg[slc] = time - time[0] 
    return reg

def make_regressor_trial_type(taxis, events, trial_type):
    # takes a dataframe of events and a time axis and returns a binary series
    reg = pd.Series(0, index=taxis, name='trial_type' + str(trial_type))
    evred = events[events['trial_type'] == trial_type]  # restrict trials
    starts = evred['trial_start']
    stops = evred['trial_over']
    pairs = zip(starts, stops)
    for p in pairs:
        reg[slice(*p)] = 1
    return reg

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

# get unique trial types
unique_trial_types = np.unique(evt['trial_type'])

# make a regressor call for each trial type
# the j=j is needed because of Python's late binding
# also, can't use every trial type (else confounded with is_in_trial, so 
# choose first trial type as baseline)
ttype_regressors = [lambda x, y, j=j: make_regressor_trial_type(x, y, j) for j in unique_trial_types[1:]]

regressor_list = [make_regressor_is_in_trial, make_regressor_is_inflating, make_regressor_elapsed_time] + ttype_regressors

regressor_frame = pd.concat([f(spks.index, evt) for f in regressor_list], 
    axis=1)