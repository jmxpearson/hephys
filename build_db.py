import MySQLdb
import pandas as pd
import pandas.io.sql as pdsql
import scipy.io as sio
import numpy as np
import h5py
import scipy.signal as ssig
import warnings

def decimate(x, decfrac, axis=-1):
    """
    just like in the scipy source code, except I use filtfilt
    and opt for hardcode defaults
    q is the fraction to decimate by; length of returned data is len(x)/q
    """
    if decfrac > 15:
        wrnstr = """You are attempting to decimate by a factor > 15. You are risking numerical instability. Consider performing multiple successive decimations instead."""
        warnings.warn(wrnstr)

    n = 8
    b, a = ssig.filter_design.cheby1(n, 0.05, 0.8/decfrac)

    y = ssig.filtfilt(b, a, x, axis=axis)

    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, decfrac)
    return y[sl]


def WriteToDB(dbname, tblname, df):
    db = MySQLdb.connect(host='localhost',
        user='root', passwd='', db=dbname)
    pdsql.write_frame(df, con=db, name=tblname,
        if_exists='append', flavor='mysql')

def SetupDB():
    # connect to database server
    db = MySQLdb.connect(host='localhost',
        user='root', passwd='')

    # create cursor
    cur = db.cursor()

    # create database and switch to it
    setupstr = """
        CREATE DATABASE bartc;
        USE bartc;
        CREATE TABLE spikes (patient TINYINT, dataset TINYINT, 
            channel TINYINT, unit TINYINT, time DECIMAL(7,3));
        CREATE TABLE lfp (patient TINYINT, dataset TINYINT, 
            channel TINYINT, unit TINYINT, time DECIMAL(7,3), 
            voltage DOUBLE);
        """
    cur.execute(setupstr)

def ImportSpikes(ftup, datadir):
    pdir = 'patient' + str(ftup[0]).zfill(3)
    fname = ('times_' + str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + 
        str(ftup[2]) + '.mat')
    fullname = datadir + pdir + '/' + fname
    dat = sio.loadmat(fullname)['cluster_class']
    
    unit = dat[:,0].astype('int')
    times = np.around(dat[:,1]/1000., decimals=3) #times are originally in ms
    sortord = np.argsort(times) #spikes aren't always time-sorted
    times = times[sortord]
    unit = unit[sortord]
    
    # now restrict to valid units:
    valid = unit == ftup[3]
    times = times[valid]
    unit = unit[valid]

    ddict = {'patient': ftup[0], 'dataset': ftup[1], 
    'channel': ftup[2], 'unit': unit, 'time': times}
    df = pd.DataFrame(ddict)

    WriteToDB('bartc', 'spikes', df)


def ImportLFP(ftup, datadir):
    pdir = 'patient' + str(ftup[0]).zfill(3)
    fname = str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + str(ftup[2]) + '.mat'
    fullname = datadir + pdir + '/' + fname
    dset = h5py.File(fullname, 'r')
    dat = dset['data'].value.squeeze()
    sr = dset['srlfp'].value.squeeze()

    vv = decimate(dat, 5)  # decimate data to 200 Hz
    sr = sr / 5;
    dt = (1. / sr).round(3)

    times = (np.arange(0, vv.shape[0]) * dt).round(3).squeeze()
    ddict = {'patient': ftup[0], 'dataset': ftup[1], 
    'channel': ftup[2], 'time': times, 'voltage': vv}

    df = pd.DataFrame(ddict)
 
    WriteToDB('bartc', 'lfp', df)

# read data
# dat = pdsql.read_frame(qstr, db)

# # write data
# pdsql.write_frame(dat, con=db, name='supp',
#     if_exists='replace', flavor='mysql')

########### tests ##################
# patient = 18
# dset = 1
# chan = 1
# unit = 2

# ftup = (patient, dset, chan, unit)


if __name__ == '__main__':

    # build database and tables
    # SetupDB()

    # locations of relevant files
    ddir = '/home/jmp33/data/bartc/plexdata/'
    spkfile = '/home/jmp33/code/hephys/valid_units.csv'

    # get list of tuples with valid channels
    ulist = []

    with open(spkfile) as infile:
        for line in infile:
            ulist.append(tuple(map(int, line.split(','))))

    # iterate through files, loading data
    # for ftup in ulist:
    #     ImportSpikes(ftup, ddir)

    ftup = 18,1,32
    ImportLFP(ftup, ddir)

