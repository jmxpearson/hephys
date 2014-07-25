from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import physutils
import hephys.dbio as dbio
import os

def make_time_frequency_plot(dtup, event_name, Tpre, Tpost, baseline_interval):

    # get lfp data
    print "Fetching data: " + str(dtup)
    lfp = dbio.fetch_all_such_LFP(dbname, *dtup).censor()

    # get events
    evt = dbio.fetch(dbname, 'events', *dtup[:2])
    times = evt[event_name].dropna()

    wav_normed, fig = lfp.avg_time_frequency(dtup[2], times, Tpre, Tpost, method='wav', normfun=physutils.norm_by_trial(baseline_interval))

    return fig

if __name__ == '__main__':

    # open data file
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

    # first, get a list of lfp channels
    setlist = pd.read_hdf(dbname, '/meta/lfplist')

    # get ready to write to file
    fname = 'pop_time_freqs.pdf'
    event_name = 'popped'
    Tpre = -1.5
    Tpost = 0.5
    baseline_interval = (-1.5, -1.35)

    # open pdf for plotting
    with PdfPages(os.path.expanduser('~/Dropbox/hephys/media/figs/' + fname)) as pdf:

        for idx, channel_inds in setlist.iterrows():
            dtup = tuple(channel_inds)
            make_time_frequency_plot(dtup, event_name, Tpre, Tpost, baseline_interval)
            titlestr = "Channel: " + str(dtup) + "\nAlign: " + event_name
            plt.title(titlestr)
            pdf.savefig()
            plt.close()
