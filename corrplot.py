"""
gist to plot cross-channel LFP correlations as a function of time 
"""
from physutils import *
from physclasses import * 
import numpy as np
import matplotlib.pyplot as plt

# pick a dataset
dtup = 17, 2 

# open data file
dbname = '/home/jmp33/data/bartc/plexdata/bartc.hdf5'

# get lfp data, normalize, downsample
df = fetch_all_such_LFP(dbname, *dtup)
df = df.zscore()
df = df.decimate(5)
    
# first, try to just look at overall correlations
# if the similarity matrix is correlation, the corresponding distance 
# measure is indeed mean-squared distance 
from scipy.spatial.distance import pdist, squareform
# have to transpose because pdist expects each ROW an observation:
Dmat = pdist(df.values.T)  
import scipy.cluster.hierarchy as clust
dendro = clust.linkage(Dmat)

# plot it:
tree = clust.dendrogram(dendro)
plt.show()

# multidimensional scaling
from sklearn import manifold
seed = 12345
mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
pos = mds.fit_transform(squareform(Dmat))
plt.scatter(pos[:, 0], pos[:, 1])
plt.show()