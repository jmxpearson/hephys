import h5py
import pandas as pd
import scipy.io as sio
import numpy as np
from physutils import decimate, MakePath

def WriteToDB(dbname, tblname, df, **kwargs):
    df.to_hdf(dbname, tblname, append=True)

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

    target = 'spikes/' + MakePath(*ftup)
    WriteToDB(datadir + 'bartc.hdf5', target, df)

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
 
    target = 'lfp/' + MakePath(*ftup) 
    WriteToDB(datadir + 'bartc.hdf5', target, df)

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
     
        target = 'censor/' + MakePath(*ftup) 
        WriteToDB(datadir + 'bartc.hdf5', target, df)

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
    isFHC = (evt[0].size < 10)

    # make a dataframe from matlab behavior, pull out event codes
    bf = pd.DataFrame(matevt)
    bf.index.names = ['trial']
    evf = bf[['ev', 'evt']]
    bf = bf.drop(['ev', 'evt'], axis=1)
    # make sure to get data types right
    bf = bf.convert_objects()
    bf['result'] = bf['result'].apply(str)
    absstart = bf['trial_start_time'][0]

    # turn each trial into a dataframe 
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

    # concatenate trials
    events = pd.concat(tlist)
    # make event names column names
    events = events.set_index('event', append=True).unstack()
    # get rid of multi-index labeling
    events.columns = pd.Index([e[1] for e in events.columns])

    # now merge both event datasets 
    df = pd.concat([bf, events], axis=1)
    df['patient'] = ftup[0]
    df['dataset'] = ftup[1]

    # if we have Plexon events, use them
    if not isFHC:
        # trial start -- sometimes a spurious event marks recording onset
        if df.shape[0] == evt[0].shape[0]:  # same number of trial starts 
            df['trial_start'] = evt[0].round(3).squeeze()
        else:
            df['trial_start'] = evt[0][1:].round(3).squeeze()

        # trial stop -- when last trial aborted, may not be present
        if df.shape[0] == evt[7].shape[0]:  # same number of trial starts 
            df['trial_over'] = evt[7].round(3).squeeze()
        else:
            df['trial_over'][:-1] = evt[7].round(3).squeeze()

        # rest of vars
        vlist = zip(['start inflating', 'stop inflating', 'banked', 'outcome', 
            'popped'], [1, 2, 3, 5, 4])
        for var in vlist:
            valid = pd.notnull(df[var[0]])
            df[valid][var[0]] = evt[var[1]].round(3).squeeze()

    # lastly, if we are missing control columns, make sure to add them
    # (important for getting schema correct on initial write)
    if not 'is_control' in df.keys():
        df['is_control'] = np.nan
    else:
        df['is_control'] = df['is_control'].astype('float64')
    if not 'ctrltime' in df.keys():
        df['ctrltime'] = np.nan
    
    # do some final tidying
    df = df[df['result'] != 'aborted']  # get rid of aborts
    target = 'events/' + MakePath(*ftup[:-1]) 
    WriteToDB(datadir + 'bartc.hdf5', target, df)


if __name__ == '__main__':

    # locations of relevant files
    ddir = '/home/jmp33/data/bartc/plexdata/'
    bdir = '/home/jmp33/data/bartc/behavior/'
    spkfile = '/home/jmp33/code/hephys/valid_units.csv'
    lfpfile = '/home/jmp33/code/hephys/lfp_channel_file.csv'
    chanfile = '/home/jmp33/code/hephys/valid_channels.csv'
    behfile = '/home/jmp33/code/hephys/behavior_file_map.csv'

    ############### events #################        
    # get all (patient, dataset) tuples from already loaded spikes
    tuplist = []
    print 'Loading Events....'
    with open(behfile) as infile:
        for line in infile:
            thistup = tuple(line.rstrip().lower().split(','))
            tuplist.append((int(thistup[0]), int(thistup[1]), thistup[2]))

    plist = [t[:-1] for t in tuplist]
    WriteToDB(ddir + 'bartc.hdf5', 'meta/evlist', 
        pd.DataFrame(plist, columns=['patient', 'dataset']))

    for ftup in tuplist:
        print ftup
        ImportEvents(ftup, ddir, bdir)

    ############### spikes #################
    # get list of tuples with valid channels
    tuplist = []
    print 'Loading Spikes....'
    with open(spkfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))

    WriteToDB(ddir + 'bartc.hdf5', 'meta/spklist', 
        pd.DataFrame(tuplist, 
            columns=['patient', 'dataset', 'channel', 'unit']))

    # iterate through files, loading data
    for ftup in tuplist:
        print ftup
        ImportSpikes(ftup, ddir)

    ############### lfp #################
    # load lfp data
    tuplist = []
    print 'Loading LFP....'
    with open(lfpfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))

    WriteToDB(ddir + 'bartc.hdf5', 'meta/lfplist', 
        pd.DataFrame(tuplist, columns=['patient', 'dataset', 'channel']))

    # iterate through files, loading data
    for ftup in tuplist:
        print ftup
        ImportLFP(ftup, ddir)

    ############### censoring #################
    # load censoring data
    tuplist = []
    print 'Loading Censoring data....'
    with open(chanfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))
    with open(lfpfile) as infile:
        for line in infile:
            tuplist.append(tuple(map(int, line.split(','))))

    WriteToDB(ddir + 'bartc.hdf5', 'meta/censlist', 
        pd.DataFrame(tuplist, columns=['patient', 'dataset', 'channel']))

    # iterate through files, loading data
    for ftup in tuplist:
        print ftup
        ImportCensor(ftup, ddir)