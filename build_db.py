import h5py
import pandas as pd
import scipy.io as sio
import numpy as np
from physutils import make_path, decimate

class DataSets:
    def __init__(self, data_dir, behavior_dir, channel_file, spk_file, lfp_file, behavior_file, output_file, plexon_event_codes, flatten_events=False):
        self.datadir = data_dir  # directory for plexon data files
        self.behdir = behavior_dir  # directory for behavior files
        self.chanfile = channel_file  # file listing spike channels to import
        self.spkfile = spk_file  # file listing spike units to import
        self.lfpfile = lfp_file  # file listing lfp channels to import
        self.behfile = behavior_file  # file mapping phys to behavior files
        self.outfile = output_file  # name of output data file
        self.plx_codes = plexon_event_codes  # dict mapping event names to plexon codes
        self.flatten_events = flatten_events  # encode events in a single column (False => 1 column per event)

    def write_to_db(self, tblname, df, **kwargs):
        df.to_hdf(self.outfile, tblname, append=True)

    def add_metadata(self, dbname, tblname, **kwargs):
        fobj = h5py.File(dbname, 'a')
        for k in kwargs:
            fobj[tblname].attrs[k] = kwargs[k]
        fobj.close()

    def import_spikes(self, ftup):
        pdir = 'patient' + str(ftup[0]).zfill(3)
        fname = ('times_' + str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + 
            str(ftup[2]) + '.mat')
        fullname = self.datadir + pdir + '/' + fname
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

        target = 'spikes/' + make_path(*ftup)
        self.write_to_db(target, df)

    def import_lfp(self, ftup):
        pdir = 'patient' + str(ftup[0]).zfill(3)
        fname = str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + str(ftup[2]) + '.mat'
        fullname = self.datadir + pdir + '/' + fname
        dset = h5py.File(fullname, 'r')
        dat = dset['data'].value.squeeze()
        sr = dset['srlfp'].value.squeeze()

        vv = dat
        decfrac = 5.0
        vv = decimate(dat, decfrac)  # decimate data to 200 Hz
        sr = sr / decfrac;
        dt = (1. / sr).round(3)

        times = (np.arange(0, vv.size) * dt).round(3).squeeze()
        ddict = {'patient': ftup[0], 'dataset': ftup[1], 
        'channel': ftup[2], 'time': times, 'voltage': vv}

        df = pd.DataFrame(ddict)
     
        target = 'lfp/' + make_path(*ftup) 
        self.write_to_db(target, df)

    def import_censor(self, ftup):
        pdir = 'patient' + str(ftup[0]).zfill(3)
        fname = (str(ftup[0]) + '.' + str(ftup[1]) + '.plx' + str(ftup[2]) + 
            '_censoring.mat')
        fullname = self.datadir + pdir + '/' + fname
        excludes = sio.loadmat(fullname)['excludes'].round(3)

        if excludes.size != 0:
            #do this in case some exclude ranges make no sense
            badrng = np.where(np.diff(excludes, axis=1) < 0)
            excludes = np.delete(excludes, badrng, axis=0)

            # get data ready
            ddict = {'patient': ftup[0], 'dataset': ftup[1], 
            'channel': ftup[2], 'start': excludes[:,0], 'stop': excludes[:,1]}

            df = pd.DataFrame(ddict)
         
            target = 'censor/' + make_path(*ftup) 
            self.write_to_db(target, df)

    def _grab_matlab_events(self, ftup):
        pdir = 'patient' + str(ftup[0]).zfill(3)
        behname = self.behdir + pdir + '/' + ftup[2]

        # load matlab events
        matevt = sio.loadmat(behname, squeeze_me=True)['data']

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
        events['event'] = events['event'].apply(str)

        return (bf, events)

    def _grab_plexon_events(self, ftup):
        pdir = 'patient' + str(ftup[0]).zfill(3)
        fname = str(ftup[0]) + '.' + str(ftup[1]) + '.plx_events.mat'
        fullname = self.datadir + pdir + '/' + fname

        # load plexon events
        evt = sio.loadmat(fullname)['evt'].squeeze()

        return evt

    def import_events(self, ftup):
        trial_variables, events = self._grab_matlab_events(ftup)
        evt = self._grab_plexon_events(ftup)

        # get number of events
        numtrials = max(events.index) + 1

        # was this an FHC recording? if so, there are no Plexon stamps
        # all events dumped into first or second slot, so other slots 
        # should have few timestamps
        isFHC = (evt[2].size < 10)

        if isFHC:  # line up events with phys
            # for now, we kludge this by just setting the clocks to be equal
            # at task start and not worrying about drift
            num_events = map(len, evt)
            startcode = np.argmax(num_events)  # these should be trial starts

            # get time of first FHC event
            FHC_start = evt[startcode][0].round(3).squeeze()

            # compensate for offset
            all_events = events.stack()
            all_events.sort()
            ephys_offset = (FHC_start - all_events.values[0]).round(3)
            events['time'] += ephys_offset

        else:  # if we have Plexon events, use them
            startcode = self.plx_codes['trial_start']
            stopcode = self.plx_codes['trial_over']

            # trial start -- sometimes a spurious event marks recording onset
            if evt[startcode].shape[0] != numtrials: 
                evt[startcode] = evt[startcode][1:]

            # trial stop -- when last trial aborted, may not be present
            if evt[stopcode].shape[0] != numtrials: 
                evt[stopcode] = np.append(evt[stopcode], np.nan)

            for var in self.plx_codes:
                # valid = pd.notnull(events[var])
                # events[var][valid] = evt[self.plx_codes[var]].round(3).squeeze()
                this_selection = events['event'] == var
                events['time'][this_selection] = evt[self.plx_codes[var]].round(3).squeeze()

        # try to make events columns: this may fail in case some events
        # can happen multiple times per trial; in that case, make each 
        # event a row and perform a join
        if self.flatten_events:
            # now merge task variables and events 
            df = events.join(trial_variables)
        else: 
            # make event names column names
            events = events.set_index('event', append=True).unstack()
            # get rid of multi-index labeling
            events.columns = pd.Index([e[1] for e in events.columns])
            # now merge task variables and events 
            df = pd.concat([trial_variables, events], axis=1)

        df['patient'] = ftup[0]
        df['dataset'] = ftup[1]

        # do some final tidying
        df = df[df['result'] != 'aborted']  # get rid of aborts
        target = 'events/' + make_path(*ftup[:-1]) 
        self.write_to_db(target, df)

    def load_all_events(self):
        # get all (patient, dataset) tuples from already loaded spikes
        tuplist = []
        print 'Loading Events....'
        with open(self.behfile) as infile:
            for line in infile:
                thistup = tuple(line.rstrip().lower().split(','))
                tuplist.append((int(thistup[0]), int(thistup[1]), thistup[2]))

        plist = [t[:-1] for t in tuplist]
        self.write_to_db('meta/evlist', 
            pd.DataFrame(plist, columns=['patient', 'dataset']))

        for ftup in tuplist:
            print ftup
            self.import_events(ftup)

    def load_all_censoring(self):
        # load censoring data
        tuplist = []
        print 'Loading Censoring data....'
        with open(self.chanfile) as infile:
            for line in infile:
                tuplist.append(tuple(map(int, line.split(','))))
        with open(self.lfpfile) as infile:
            for line in infile:
                tuplist.append(tuple(map(int, line.split(','))))

        self.write_to_db('meta/censlist', 
            pd.DataFrame(tuplist, columns=['patient', 'dataset', 'channel']))

        # iterate through files, loading data
        for ftup in tuplist:
            print ftup
            self.import_censor(ftup)

    def load_all_spikes(self):
        # get list of tuples with valid channels
        tuplist = []
        print 'Loading Spikes....'
        with open(self.spkfile) as infile:
            for line in infile:
                tuplist.append(tuple(map(int, line.split(','))))

        self.write_to_db('meta/spklist', 
            pd.DataFrame(tuplist, 
                columns=['patient', 'dataset', 'channel', 'unit']))

        # iterate through files, loading data
        for ftup in tuplist:
            print ftup
            self.import_spikes(ftup)

    def load_all_lfp(self):
        # load lfp data
        tuplist = []
        print 'Loading LFP....'
        with open(self.lfpfile) as infile:
            for line in infile:
                tuplist.append(tuple(map(int, line.split(','))))

        self.write_to_db('meta/lfplist', 
            pd.DataFrame(tuplist, columns=['patient', 'dataset', 'channel']))

        # iterate through files, loading data
        for ftup in tuplist:
            print ftup
            self.import_lfp(ftup)

    def load_all(self):
        self.load_all_events()
        self.load_all_censoring()
        self.load_all_spikes()
        self.load_all_lfp()



if __name__ == '__main__':

    # locations of relevant files
    ddir = '/home/jmp33/data/bartc/plexdata/'
    bdir = '/home/jmp33/data/bartc/behavior/'
    chanfile = '/home/jmp33/code/hephys/bartc/valid_channels.csv'
    spkfile = '/home/jmp33/code/hephys/bartc/valid_units.csv'
    lfpfile = '/home/jmp33/code/hephys/bartc/lfp_channel_file.csv'
    behfile = '/home/jmp33/code/hephys/bartc/behavior_file_map.csv'
    outfile = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

    plexon_event_codes = {'start inflating': 1, 'stop inflating': 2, 'banked': 3, 'outcome': 5, 'popped': 4, 'trial_start': 0, 'trial_over': 7}

    dset = DataSets(ddir, bdir, chanfile, spkfile, lfpfile, behfile, outfile, plexon_event_codes)

    dset.load_all()