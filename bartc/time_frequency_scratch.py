from physutils import *
from physclasses import *
import os

os.chdir(os.path.expanduser('~/code/hephys/bartc'))
dtup = 17, 2, 17

# open data file
dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

# get lfp data
print "Fetching Data..."
lfp = fetch_all_such_LFP(dbname, *dtup).censor()

# get events
evt = fetch(dbname, 'events', *dtup[:2])
stops = evt['stop inflating'].dropna()

Tpre = -1.5
Tpost = 0.5
window = 0.250
percent_overlap = .95

spec = lfp.avg_time_frequency(17, stops, Tpre, Tpost, method='spec', winlen=window, frac_overlap=percent_overlap)

wav = lfp.avg_time_frequency(17, stops, Tpre, Tpost, method='wav')

norm_range = (-1.5, -1.35)
wav_normed, im = lfp.avg_time_frequency(17, stops, Tpre, Tpost, method='wav', normfun=norm_by_trial(norm_range))