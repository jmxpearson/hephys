"""
gist to plot cross-channel LFP correlations as a function of time 
"""
from physutils import *
from physclasses import * 
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as clust
from sklearn import manifold
import os

# pick a dataset
dtup = 23, 1 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# get lfp data, normalize, downsample
print 'Fetching data...'
df = fetch_all_such_LFP(dbname, *dtup)
df = df.zscore()
df = df.decimate(5)
    
# first, try to just look at overall correlations
# if the similarity matrix is correlation, the corresponding distance 
# measure is indeed mean-squared distance 
# have to transpose because pdist expects each ROW an observation:
print 'Calculating dendrogram...'
Dmat = pdist(df.values.T)  
dendro = clust.linkage(Dmat)

# plot it:
tree = clust.dendrogram(dendro)
plt.show()

# multidimensional scaling
seed = 12345
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit_transform(squareform(Dmat))
plt.scatter(pos[:, 0], pos[:, 1])
plt.show()

############################ make a movie! ###############
# set some parameters
Tchunk = 5  # size of moving window (in s)
overlap = 0.95  # percent overlap
Nchunk = np.ceil(Tchunk * df.meta['sr']).astype('int')
Nshift = np.ceil((1 - overlap) * Nchunk).astype('int')

# set up initial MDS 
seed = 12345
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity='precomputed', n_jobs=1)

# set up plotting:
files = []
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
outdir = '/home/jmp33/data/bartc/movies'

offset = 0
counter = 0
while (offset + Nchunk < len(df.index)):
    if counter % 100 == 0:
        print 'counter = ' + str(counter)
    chunk = df.iloc[offset:offset + Nchunk]
    dist = squareform(pdist(chunk.values.T))
    if offset == 0:
        pos = mds.fit_transform(squareform(pdist(chunk.values.T)))
        # coords = [pos]

        # now reset parameters of mds such that in future, we *have to* start
        # from inits (and use only that start)
        mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity='precomputed', n_jobs=1, n_init=1)
    else:
        pos = mds.fit_transform(squareform(pdist(chunk.values.T)), init=oldpos)
        # coords.append(pos)

    oldpos = pos
    offset += Nshift
    counter += 1

    # plot and save image to file
    ax.cla()
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    ax.scatter(pos[:, 0], pos[:, 1])
    fname = os.path.join(outdir, '_tmp' + str(counter) + '.png')
    fig.savefig(fname)
    files.append(fname)

print 'Making movie animation.mpg - this make take a while'
curdir = os.getcwd() 
os.chdir(outdir)
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o " + '.'.join(map(str, dtup)) + ".mpg")
os.system("rm _tmp*")
os.chdir(curdir)