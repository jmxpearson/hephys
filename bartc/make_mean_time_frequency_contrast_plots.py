"""
For each dataset, make a time-frequency plot that is a 
contrast between similar plots for two conditions. The plot for 
each condition is a mean across channels within the dataset.
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import physutils
import hephys.dbio as dbio
import os

def make_time_frequency_plot(dtup, event_names, Tpre, Tpost, freqs, baseline_interval, thresh):

    # get lfp data
    print "Fetching data: " + str(dtup)
    lfp = dbio.fetch_all_such_LFP(dbname, *dtup).censor().zscore()

    # get events
    evt = dbio.fetch(dbname, 'events', *dtup[:2])
    times0 = evt[event_names[0]].dropna()
    times1 = evt[event_names[1]].dropna()

    # horrible kludge to drop pathological channels
    bad_channel_list = [(16, 2, 22), (18, 1, 32)]

    all_wavs0 = []
    all_wavs1 = []
    for channel in lfp.columns:
        # horrible kludge to exclude pathological channels
        if dtup + (channel,) in bad_channel_list: 
            continue
        print "Channel " + str(channel)
        wav_normed0, im = lfp.avg_time_frequency(channel, times0, Tpre, Tpost, method='wav', doplot=False, normfun=None, freqs=freqs)
        all_wavs0.append(wav_normed0)
        wav_normed1, im = lfp.avg_time_frequency(channel, times1, Tpre, Tpost, method='wav', doplot=False, normfun=None, freqs=freqs)
        all_wavs1.append(wav_normed1)

    # take mean power across all channels
    all_wav_mean0 = reduce(lambda x, y: x.add(y, fill_value=0), all_wavs0) / len(all_wavs0)
    all_wav_mean1 = reduce(lambda x, y: x.add(y, fill_value=0), all_wavs1) / len(all_wavs1)

    fig1 = physutils.tf.plot_time_frequency(all_wav_mean0/all_wav_mean1)
    return (fig1,) 

if __name__ == '__main__':

    # open data file
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

    # first, get a list of lfp channels
    setlist = pd.read_hdf(dbname, '/meta/lfplist')[['patient', 'dataset']].drop_duplicates()

    # get ready to write to file
    # fname = 'stop_vs_start.pdf'
    fname = 'stop_vs_pop.pdf'
    # event_names = ['stop inflating', 'start inflating']
    event_names = ['stop inflating', 'popped']
    Tpre = -1.5
    Tpost = 0.5
    baseline_interval = (-1.5, -1.35)
    freqs = np.exp(np.linspace(np.log(2.5), np.log(50)))
    thresh = 1.5

    # open pdf for plotting
    with PdfPages(os.path.expanduser('~/Dropbox/hephys/media/figs/' + fname)) as pdf:

        for _, channel_inds in setlist.iterrows():
            dtup = tuple(channel_inds)
            figs = make_time_frequency_plot(dtup, event_names, Tpre, Tpost, freqs, baseline_interval, thresh)
            titlestr = "Channel: " + str(dtup) + "\nAlign: " + event_names[0] + " / " + event_names[1]
            for f in figs:
                f.suptitle(titlestr)
                pdf.savefig(f)
                plt.close(f)
