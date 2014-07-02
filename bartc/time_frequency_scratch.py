from physutils import *
from physclasses import *
import os

os.chdir(os.path.expanduser('~/code/hephys/bartc'))
dtup = 17, 2

# open data file
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# get lfp data
print "Fetching Data..."
lfp = fetch_all_such_LFP(dbname, *dtup)

# get events
evt = fetch(dbname, 'events', *dtup)
stops = evt['stop inflating'].dropna()

# brabe a particular channel
channel = lfp.dataframe[17]

Tpre = -1.5
Tpost = 0.5
window = 0.250
percent_overlap = .95

spec = avg_spectrogram(channel, window, percent_overlap, stops, Tpre, Tpost)

plot_spectrogram(spec)