from build_db import *

# locations of relevant files
ddir = '/home/jmp33/data/gonogo/plexdata/'
bdir = '/home/jmp33/data/gonogo/behavior/'
chanfile = '/home/jmp33/code/hephys/gonogo/valid_channels.csv'
spkfile = '/home/jmp33/code/hephys/gonogo/valid_units.csv'
lfpfile = '/home/jmp33/code/hephys/gonogo/lfp_channel_file.csv'
behfile = '/home/jmp33/code/hephys/gonogo/behavior_file_map.csv'
outfile = '/home/jmp33/data/gonogo/plexdata/gonogo.hdf5'

# plexon_event_codes = {'start inflating': 1, 'stop inflating': 2, 'banked': 3, 'outcome': 5, 'popped': 4, 'trial_start': 0, 'trial_over': 7}
plexon_event_codes = {}
flatten_events = True

dset = DataSets(ddir, bdir, chanfile, spkfile, lfpfile, behfile, outfile, plexon_event_codes, flatten_events)

dset.load_all()