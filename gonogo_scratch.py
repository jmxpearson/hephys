from physutils import *
from pylab import *

def plot_psth(dbname, ftup, evt_name):
    df = fetch(dbname, 'spikes', *ftup) 
    evt = fetch(dbname, 'events', *ftup[:2])

    binsize = 0.050  # 50 ms bin
    binned = binspikes(df, binsize)

    evt_set = evt[evt['event'] == evt_name]['time']

    psth = evtsplit(binned, evt_set, -0.5, 1).mean(axis=1) / binsize

    smpsth = smooth(psth, 0.4)

    plot(smpsth.index, smpsth)
    show()

if __name__ == '__main__':
    dbname = '/home/jmp33/data/gonogo/plexdata/gonogo.hdf5'
    patient = 33 
    dataset = 2
    channel = 1
    unit = 1
    evt_name = 'no_response'

    plot_psth(dbname, (patient, dataset, channel, unit), evt_name)