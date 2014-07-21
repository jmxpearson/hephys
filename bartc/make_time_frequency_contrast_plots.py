from matplotlib.backends.backend_pdf import PdfPages
from physclasses import *
from physutils import *
import os

def make_time_frequency_plot(dtup, event_names, Tpre, Tpost, baseline_interval):

    # get lfp data
    print "Fetching data: " + str(dtup)
    lfp = fetch_all_such_LFP(dbname, *dtup).censor()

    # get events
    evt = fetch(dbname, 'events', *dtup[:2])
    times0 = evt[event_names[0]].dropna()
    times1 = evt[event_names[1]].dropna()

    wav_normed0, fig = lfp.avg_time_frequency(dtup[2], times0, Tpre, Tpost, method='wav', doplot=False, normfun=norm_by_trial(baseline_interval))
    wav_normed1, fig = lfp.avg_time_frequency(dtup[2], times1, Tpre, Tpost, method='wav', doplot=False, normfun=norm_by_trial(baseline_interval))

    fig = plot_time_frequency(wav_normed0 / wav_normed1) 

    return fig

if __name__ == '__main__':

    os.chdir(os.path.expanduser('~/code/hephys/bartc'))

    # open data file
    dbname = os.path.expanduser('~/data/bartc/plexdata/bartc.hdf5')

    # first, get a list of lfp channels
    setlist = pd.read_hdf(dbname, '/meta/lfplist')

    # get ready to write to file
    fname = 'stop_vs_pop_time_freqs.pdf'
    event_names = ['popped', 'stop inflating']
    Tpre = -1.5
    Tpost = 0.5
    baseline_interval = (-1.5, -1.35)

    # open pdf for plotting
    with PdfPages(os.path.expanduser('~/Dropbox/hephys/media/figs/' + fname)) as pdf:

        for idx, channel_inds in setlist.iterrows():
            dtup = tuple(channel_inds)
            make_time_frequency_plot(dtup, event_names, Tpre, Tpost, baseline_interval)
            titlestr = "Channel: " + str(dtup) + "\nAlign: " + event_names[0] + " / " + event_names[1]
            plt.title(titlestr)
            plt.ylim((0, 40))
            pdf.savefig()
            plt.close()
