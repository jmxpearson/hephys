from physutils import *
from physclasses import *
import os

os.chdir(os.path.expanduser('~/code/hephys/bartc'))
dtup = 17, 2, 17

# open data file
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# get lfp data
print "Fetching Data..."
lfp = fetch_all_such_LFP(dbname, *dtup)

# get events
evt = fetch(dbname, 'events', *dtup[:2])
stops = evt['stop inflating'].dropna()

# brabe a particular channel
channel = lfp

Tpre = -1.5
Tpost = 0.5
window = 0.250
percent_overlap = .95

spec = avg_time_frequency(channel, spectrogram, stops, Tpre, Tpost, window, percent_overlap) 

plot_time_frequency(spec)

tf = avg_time_frequency(channel, continuous_wavelet, stops, Tpre, Tpost)
tfint = interpolate_time_frequency(tf)

plot_time_frequency(tfint)

