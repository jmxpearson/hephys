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
    lfp = dbio.fetch_all_such_LFP(dbname, *dtup).censor()

    # get events
    evt = dbio.fetch(dbname, 'events', *dtup[:2])
    times0 = evt[event_names[0]].dropna()
    times1 = evt[event_names[1]].dropna()

    nf = physutils.norm_by_trial(baseline_interval, method='subtraction')

    contr_tf, fig1 = lfp.contrast_time_frequency(dtup[2], [times0, times1], Tpre, Tpost, method='wav', normfun=nf, doplot=True, freqs=freqs)

    mcontr, fig2 = lfp.significant_time_frequency(dtup[2], [times0, times1], Tpre, Tpost, thresh=thresh, niter=1000, method='wav', doplot=True, normfun=nf, freqs=freqs)

    return fig1, fig2

if __name__ == '__main__':

    # open data file
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

    # first, get a list of lfp channels
    setlist = pd.read_hdf(dbname, '/meta/lfplist')

    # get ready to write to file
    fname = 'stop_vs_pop_time_freqs.pdf'
    event_names = ['stop inflating', 'popped']
    Tpre = -1.5
    Tpost = 0.5
    baseline_interval = (-1.5, -1.35)
    freqs = np.exp(np.linspace(np.log(2.5), np.log(50)))
    thresh = 1.

    # open pdf for plotting
    with PdfPages(os.path.expanduser('~/Dropbox/hephys/media/figs/' + fname)) as pdf:

        for idx, channel_inds in setlist.iterrows():
            dtup = tuple(channel_inds)
            figs = make_time_frequency_plot(dtup, event_names, Tpre, Tpost, freqs, baseline_interval, thresh)
            titlestr = "Channel: " + str(dtup) + "\nAlign: " + event_names[0] + " / " + event_names[1]
            for f in figs:
                f.suptitle(titlestr)
                # f.gca().set_ylim(0, 40)
                pdf.savefig(f)
                plt.close(f)
