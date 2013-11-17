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
        CREATE TABLE censor (patient TINYINT, dataset TINYINT, 
            channel TINYINT, start DECIMAL(7,3), 
            stop DECIMAL(7,3));
        CREATE TABLE events (patient TINYINT, dataset TINYINT, 
            trial_start_time DECIMAL(7,3), this_balloon TINYINT,
            trial_type TINYINT, points INT, inflate_time DECIMAL(7,3),
            this_run TINYINT, rt DECIMAL(7,3), score INT, 
            banked DECIMAL(7,3), outcome DECIMAL(7,3), popped DECIMAL(7,3),
            start_inflating DECIMAL(7,3), stop_inflating DECIMAL(7,3),
            trial_over DECIMAL(7,3), trial_start DECIMAL(7,3));
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

    times = (np.arange(0, vv.size) * dt).round(3).squeeze()
    ddict = {'patient': ftup[0], 'dataset': ftup[1], 
    'channel': ftup[2], 'time': times, 'voltage': vv}

    df = pd.DataFrame(ddict)
 
    WriteToDB('bartc', 'lfp', df)

def ImportCensor(ftup, datadir):
    pdir = 'patient' + str(ftup[0]).zfill(3)
    fname = (str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + str(ftup[2]) + 
        '_censoring.mat')
    fullname = datadir + pdir + '/' + fname
    excludes = sio.loadmat(fullname)['excludes'].round(3)

    if excludes.size != 0:
        ddict = {'patient': ftup[0], 'dataset': ftup[1], 
        'channel': ftup[2], 'start': excludes[:,0], 'stop': excludes[:,1]}

        df = pd.DataFrame(ddict)
     
        WriteToDB('bartc', 'censor', df)

def ImportEvents(ftup, datadir, behdir):
    pdir = 'patient' + str(ftup[0]).zfill(3)
    fname = str(ftup[0]) + '.' + str(ftup[1]) + '.plx_events.mat' 
    fullname = datadir + pdir + '/' + fname
    behname = behdir + pdir + '/' + ftup[2]

    # load plexon events
    evt = sio.loadmat(fullname)['evt'].squeeze()

    # load matlab events
    matevt = sio.loadmat(behname, squeeze_me = True)['data']

    # was this an FHC recording? if so, there are no Plexon stamps
    isFHC = (evt.size < 8) or (evt[-1].size == 1)

    # make a dataframe from matlab behavior, pull out event codes
    bf = pd.DataFrame(matevt)
    bf.index.names = ['trial']
    evf = bf[['ev', 'evt']]
    bf = bf.drop(['ev', 'evt'], axis=1)
    absstart = bf['trial_start_time'][0]

    # get all unique event types
    tlist = []
    for ind in evf.index:
        thistrial = {'event': evf['ev'][ind], 'time': evf['evt'][ind]}
        miniframe = pd.DataFrame(thistrial, 
            index=pd.Index(ind * np.ones_like(evf['ev'][ind]), name='trial'))
        # now take care of timestamps by putting times in seconds and adjusting
        # by start time of first trial
        miniframe['time'] /= 1000. 
        miniframe['time'] += bf['trial_start_time'][ind] - absstart
        miniframe['time'] = miniframe['time'].round(3)  # round to ms
        tlist.append(miniframe)
    events = pd.concat(tlist)
    # make event names column names
    events = events.set_index('event', append=True).unstack()
    # get rid of multi-index labeling
    events.columns = pd.Index([e[1] for e in events.columns])

    # now merge all data for each trial
    df = pd.concat([bf, events], axis=1)
    df['patient'] = ftup[0]
    df['dataset'] = ftup[1]

    # if we have Plexon events, use them
    if not isFHC:
        # trial start
        if df.shape[0] == evt[0].shape[0]:  # same number of trial starts 
            df['trial_start'] = evt[0].round(3)
        else:
            df['trial_start'] = evt[0][1:].round(3)

        # trial stop
        if df.shape[0] == evt[7].shape[0]:  # same number of trial starts 
            df['trial_over'] = evt[7].round(3)
        else:
            df['trial_over'][:-1] = evt[7].round(3)

        # rest of vars
        vlist = zip(['start inflating', 'stop inflating', 'banked', 'outcome', 
            'popped'], [1, 2, 3, 5, 4])
        for var in vlist:
            valid = pd.notnull(df[var[0]])
            df[valid][var[0]] = evt[var[1]].round(3).squeeze()

    # do some final tidying
    df = df[df['result'] != 'aborted']  # get rid of aborts
    df = df.where((pd.notnull(df)), None)  # replace NaN with None

    WriteToDB('bartc', 'events', df)

# read data
# dat = pdsql.read_frame(qstr, db)

# # write data
# pdsql.write_frame(dat, con=db, name='supp',
#     if_exists='replace', flavor='mysql')

if __name__ == '__main__':

    # build database and tables
    SetupDB()

    # locations of relevant files
    ddir = '/home/jmp33/data/bartc/plexdata/'
    bdir = '/home/jmp33/data/bartc/behavior/'
    spkfile = '/home/jmp33/code/hephys/valid_units.csv'
    lfpfile = '/home/jmp33/code/hephys/lfp_channel_file.csv'
    chanfile = '/home/jmp33/code/hephys/valid_channels.csv'
    behfile = '/home/jmp33/code/hephys/behavior_file_map.csv'

    ############### spikes #################
    # get list of tuples with valid channels
    tuplist = []
    print 'Loading Spikes....'
    with open(spkfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))

    # iterate through files, loading data
    for ftup in tuplist:
        print ftup
        ImportSpikes(ftup, ddir)

    ############### events #################        
    # get all (patient, dataset) tuples from already loaded spikes
    tuplist = []
    print 'Loading Events....'
    with open(behfile) as infile:
        for line in infile:
            thistup = tuple(line.rstrip().lower().split(','))
            tuplist.append((int(thistup[0]), int(thistup[1]), thistup[2]))

    for ftup in tuplist:
        print ftup
        ImportEvents(ftup, ddir, bdir)

    ############### lfp #################
    # load lfp data
    tuplist = []
    print 'Loading LFP....'
    with open(lfpfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))

    # iterate through files, loading data
    for ftup in tuplist:
        print ftup
        ImportLFP(ftup, ddir)

    ############### censoring #################
    # load lfp data
    tuplist = []
    print 'Loading Censoring data....'
    with open(chanfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))

    # iterate through files, loading data
    for ftup in tuplist:
        print ftup
        ImportCensor(ftup, ddir)