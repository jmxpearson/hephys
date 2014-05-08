from physutils import *
from pylab import *

# get all spikes for a given unit 
dbname = '/home/jmp33/data/gonogo/plexdata/gonogo.hdf5'
patient = 32
dataset = 1
channel = 1
unit = 1
df = fetch(dbname, 'spikes', patient, dataset, channel, unit)
evt = fetch(dbname, 'events', patient, dataset)

binsize = 0.050  # 50 ms bin
binned = binspikes(df, binsize)

event_name = 'responded'
evt_set = evt[evt['event'] == event_name]['time']

psth = evtsplit(binned, evt_set, -1, 1).mean(axis=1)

smpsth = smooth(psth, 0.4)

plot(smpsth.index, smpsth)
show()