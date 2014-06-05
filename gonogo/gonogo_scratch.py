from physutils import *
from pylab import *

def plot_psth(dbname, ftup, evt_name):
    df = fetch(dbname, 'spikes', *ftup) 
    evt = fetch(dbname, 'events', *ftup[:2])

    binsize = 0.050  # 50 ms bin
    binned = binspikes(df, binsize)

    evt_set = evt[evt['event'] == evt_name]['time']

    psth = evtsplit(binned, evt_set, -0.5, 1).mean(axis=1) / binsize
    cvpsth = evtsplit(binned, evt_set, -0.5, 1).std(axis=1) / (binsize * sqrt(evt_set.shape[0]))

    smpsth = smooth(psth, 0.4)
    smcv = smooth(cvpsth, 0.4)

    plot(smpsth.index, smpsth, color='k')
    plot(smpsth.index, smpsth + smcv, color='k', linestyle='--')
    plot(smpsth.index, smpsth - smcv, color='k', linestyle='--')
    show()

def plot_gonogo_contrast(dbname, ftup, align):
    df = fetch(dbname, 'spikes', *ftup)
    evt = fetch(dbname, 'events', *ftup[:2])

    binsize = 0.050  # 50 ms bin
    binned = binspikes(df, binsize)

    did_go = evt['result'] == 'hit'
    evt_go = evt[did_go & (evt['event'] == align)]['time']
    evt_nogo = evt[~did_go & (evt['event'] == align)]['time']

    psth_go = evtsplit(binned, evt_go, -0.5, 1).mean(axis=1) / binsize
    psth_nogo = evtsplit(binned, evt_nogo, -0.5, 1).mean(axis=1) / binsize

    cv_go = evtsplit(binned, evt_go, -0.5, 1).std(axis=1) / (binsize * sqrt(evt_go.shape[0]))
    cv_nogo = evtsplit(binned, evt_nogo, -0.5, 1).std(axis=1) / (binsize * sqrt(evt_nogo.shape[0]))

    psth_go = smooth(psth_go, 0.4)
    psth_nogo = smooth(psth_nogo, 0.4)

    smcv_go = smooth(cv_go, 0.4)
    smcv_nogo = smooth(cv_nogo, 0.4)

    plot(psth_go.index, psth_go, color='g')
    plot(psth_go.index, psth_go + smcv_go, color='g', linestyle='--')
    plot(psth_go.index, psth_go - smcv_go, color='g', linestyle='--')
    plot(psth_nogo.index, psth_nogo, color='r')
    plot(psth_nogo.index, psth_nogo + smcv_nogo, color='r', linestyle='--')
    plot(psth_nogo.index, psth_nogo - smcv_nogo, color='r', linestyle='--')
    show()

if __name__ == '__main__':
    dbname = '/home/jmp33/data/gonogo/plexdata/gonogo.hdf5'
    patient = 33 
    dataset = 1
    channel = 1
    unit = 1
    ftup = (patient, dataset, channel, unit)
    evt_name = 'no_response'

    plot_psth(dbname, ftup, evt_name)
    plot_gonogo_contrast(dbname, ftup, 'trial_start')

