from build_db import *

# locations of relevant files
ddir = '/home/jmp33/data/gonogo/plexdata/'
bdir = '/home/jmp33/data/gonogo/behavior/'
chanfile = '/home/jmp33/code/hephys/gonogo/data/valid_channels.csv'
spkfile = '/home/jmp33/code/hephys/gonogo/data/valid_units.csv'
lfpfile = '/home/jmp33/code/hephys/gonogo/data/lfp_channel_file.csv'
behfile = '/home/jmp33/code/hephys/gonogo/data/behavior_file_map.csv'
outfile = '/home/jmp33/data/gonogo/plexdata/gonogo.hdf5'

plexon_event_codes = {}
flatten_events = True

dset = DataSets(ddir, bdir, chanfile, spkfile, lfpfile, behfile, outfile, plexon_event_codes, flatten_events)

dset.load_all()