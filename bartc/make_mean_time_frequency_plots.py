"""
For all datasets, calculate the peri-event time-frequency plot, averaged
across all channels.
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import physutils
import hephys.dbio as dbio
import os

def make_time_frequency_plot(dtup, event_name, Tpre, Tpost, freqs, baseline_interval):

    # get lfp data
    print "Fetching data: " + str(dtup)
    lfp = dbio.fetch_all_such_LFP(dbname, *dtup).censor().zscore()

    # get events
    evt = dbio.fetch(dbname, 'events', *dtup[:2])
    times = evt[event_name].dropna()

    # horrible kludge to drop pathological channels
    bad_channel_list = [(16, 2, 22), (18, 1, 32)]

    all_wavs = []
    for channel in lfp.columns:
        # horrible kludge to exclude (16, 2, 22)
        if dtup + (channel,) in bad_channel_list: 
            continue
        print "Channel " + str(channel)
        wav_normed, im = lfp.avg_time_frequency(channel, times, Tpre, Tpost, method='wav', doplot=False, normfun=None)
        all_wavs.append(wav_normed)

    # take mean power across all channels
    all_wav_mean = reduce(lambda x, y: x.add(y, fill_value=0), all_wavs) / len(all_wavs)

    # normalize across frequencies
    normfun = lambda x: x / x[slice(*baseline_interval)].mean()

    fig = physutils.plot_time_frequency(normfun(all_wav_mean))

    return fig

if __name__ == '__main__':

    os.chdir(os.path.expanduser('~/code/hephys/bartc'))

    # open data file
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

    # first, get a list of lfp datasets
    setlist = pd.read_hdf(dbname, '/meta/lfplist')[['patient', 'dataset']].drop_duplicates()

    # get ready to write to file
    fname = 'mean_time_freqs.pdf'
    event_name = 'stop inflating'
    Tpre = -1.5
    Tpost = 0.5
    baseline_interval = (-1.5, -1.35)
    freqs = np.exp(np.linspace(np.log(2.5), np.log(50)))

    # open pdf for plotting
    with PdfPages(os.path.expanduser('~/Dropbox/hephys/media/figs/' + fname)) as pdf:

        for idx, channel_inds in setlist.iterrows():
            dtup = tuple(channel_inds)
            make_time_frequency_plot(dtup, event_name, Tpre, Tpost, freqs, baseline_interval)
            titlestr = "Channel: " + str(dtup) + "\nAlign: " + event_name
            plt.title(titlestr)
            pdf.savefig()
            plt.close()
