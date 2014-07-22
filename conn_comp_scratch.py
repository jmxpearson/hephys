import os
os.chdir('/home/jmp33/code/hephys/bartc')

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import numpy as np
from physutils import *

orig_img = mpimg.imread('circles.png')

def rgb2gray(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

    img = np.around(rgb2gray(orig_img), 3)
    plt.imshow(img, cmap = cm.Greys_r)
    plt.show()

bimg = rgb2gray(orig_img) < 0.9  # treat code black as True

clust_map = label_clusters(bimg)
plt.imshow(clust_map, cmap=cm.Greys)
plt.show()

aa, idx = np.unique(clust_map, return_inverse=True)