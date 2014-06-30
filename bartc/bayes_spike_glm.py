import numpy as np
import pandas as pd
import os
import pystan

def make_data_dict(df):
    # munge dataframe into dictionary of data for stan
    data_dict = {}
    data_dict['N'] = df.shape[0]
    data_dict['P'] = df.shape[1] - 2  # don't count 'count' or index
    data_dict['y'] = df['count']
    data_dict['X'] = df.drop(['time', 'count'], axis=1)
    return data_dict

def stan_to_series(fit, data_dict):
    # fit is a stan point estimation object
    names = data_dict['X'].columns
    return pd.Series(fit['beta'], index=names)

# set a random seed
np.random.seed(12345)

# name of database to use
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')
datadir = '/home/jmp33/data/bartc/'

# first, get a list of lfp channels
setlist = pd.read_hdf(dbname, '/meta/spklist')

# compile Stan model -- only need to do this once
sm = pystan.StanModel(file='spike_glm.stan')

coeffs = pd.DataFrame()

for idx, row in setlist.iterrows():
    dtup = tuple(row) 
    datfile = datadir + '.'.join(map(str, dtup)) + '.spkglmdata.csv'
    df = pd.read_csv(datfile)

    data_dict = make_data_dict(df)

    # for now, do a simple point estimate via optimization
    fit = sm.optimizing(data=data_dict)
    coeffs = coeffs.append(stan_to_series(fit, data_dict), ignore_index=True)

perc_change = np.exp(coeffs) * 100 - 100
nonshrunk = abs(perc_change) > 10
print nonshrunk.sum()

import matplotlib.pyplot as plt
plt.subplot(1, 2, 1)
plt.imshow(perc_change, interpolation='nearest')
plt.subplot(1, 2, 2)
plt.imshow(nonshrunk, interpolation='nearest')
plt.show()