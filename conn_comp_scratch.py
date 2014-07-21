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

bimg = rgb2gray(orig_img) < 0.999  # treat code black as True

clust_map = np.zeros(bimg.shape)
uf = UnionFind()

it = np.nditer(bimg, flags=['multi_index'])
while not it.finished:
    idx = it.multi_index

    # now traverse matrix by rows
    # if present cell is not background
    if bimg[idx]:
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

# second pass: label by root
it = np.nditer(bimg, flags=['multi_index'])
while not it.finished:
    idx = it.multi_index

    if bimg[idx]:
        clust_map[idx] = uf.find(idx)[0]

    it.iternext()
print "foo"