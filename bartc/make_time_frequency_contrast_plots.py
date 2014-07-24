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

    contr_tf, fig1 = lfp.contrast_time_frequency(dtup[2], [times0, times1], Tpre, Tpost, method='wav', normfun=norm_by_trial(baseline_interval), 
        doplot=True)

    mcontr, fig2 = lfp.significant_time_frequency(dtup[2], [times0, times1], Tpre, 
        Tpost, thresh=(2.5,), niter=1000, method='wav', doplot=True, 
        normfun=norm_by_trial(baseline_interval))

    return fig1, fig2

if __name__ == '__main__':

    os.chdir(os.path.expanduser('~/code/hephys/bartc'))

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

    # open pdf for plotting
    with PdfPages(os.path.expanduser('~/Dropbox/hephys/media/figs/' + fname)) as pdf:

        for idx, channel_inds in setlist.iterrows():
            dtup = tuple(channel_inds)
            figs = make_time_frequency_plot(dtup, event_names, Tpre, Tpost, baseline_interval)
            titlestr = "Channel: " + str(dtup) + "\nAlign: " + event_names[0] + " / " + event_names[1]
            for f in figs:
                f.suptitle(titlestr)
                f.gca().set_ylim(0, 40)
                pdf.savefig(f)
                plt.close(f)
            # plt.title(titlestr)
            # plt.ylim((0, 40))
            # pdf.savefig()
            # plt.close()
