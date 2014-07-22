import os
os.chdir('/home/jmp33/code/hephys/bartc')

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import numpy as np
from unionfind import UnionFind

orig_img = mpimg.imread('circles.png')

def rgb2gray(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

    img = np.around(rgb2gray(orig_img), 3)
    plt.imshow(img, cmap = cm.Greys_r)
    plt.show()

bimg = rgb2gray(orig_img) < 0.9  # treat code black as True

def label_clusters(img):
    """
    Given a numpy array of Booleans, labels all False entries with 0 and
    all True entries within a connected component by a unique integer.
    Assumes a 4-neighborhood for connectivity. Returns labeled array
    of same shape as input.
    """
    clust_map = np.zeros(img.shape)
    uf = UnionFind()

    # first loop: traverse image by pixels, constructing union-find
    # for connected components
    it = np.nditer(img, flags=['multi_index'])
    while not it.finished:
        idx = it.multi_index

        # if present cell is not background
        if img[idx]:
            uf.add(idx)  # add to union-find
            if idx[0] > 0:
                left = (idx[0] - 1, idx[1])
                if uf.find(left):
                    uf.union(left, idx)  # attach to left neighbor
            if idx[1] > 0:
                above = (idx[0], idx[1] - 1)
                if uf.find(above):
                    uf.union(above, idx)  # attach to upper neighbor

        it.iternext()

    # get roots of union-find, construct a code dict to relabel them
    # as integers
    roots = set(map(lambda x: x[0], uf.nodes.values()))
    code_dict = dict(zip(roots, np.arange(1, len(roots))))

    # second pass: label by root
    it = np.nditer(img, flags=['multi_index'])
    while not it.finished:
        idx = it.multi_index

        if img[idx]:
            clust_map[idx] = code_dict[uf.find(idx)[0]]

        it.iternext()

    return clust_map

clust_map = label_clusters(bimg)
