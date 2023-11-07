#!/usr/bin/env python 
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
# load data 

#iter1 = np.loadtxt("iteration1_data.csv", delimiter=",", skiprows=1)
#
## convolve get temperatures 
#window = 2000
#window_avg = np.convolve(iter1[:,1], np.ones((window,))/window, mode="same")
#
#steps = iter1[:,0]
#energy_std = iter1[:,2]
#
#
##https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
#points = np.array([steps,energy_std]).T.reshape(-1,1,2)
#segments = np.concatenate([points[:-1],points[1:]],axis=1)
#fig,ax = plt.subplots()
#
#norm = plt.Normalize(2800,4000)
#lc = LineCollection(segments,cmap="magma",norm=norm)
#lc.set_array(window_avg)
#lc.set_linewidth(1)
#line = ax.add_collection(lc)
#
##ax.set_ylim([-0.001,0.5])
#plt.savefig("test.png")

x = np.linspace(0, 3 * np.pi, 500)
y = np.sin(x)
dydx = np.cos(0.5 * (x[:-1] + x[1:]))  # first derivative

# Create a set of line segments so that we can color them individually
# This creates the points as an N x 1 x 2 array so that we can stack points
# together easily to get the segments. The segments array for line collection
# needs to be (numlines) x (points per line) x 2 (for x and y)
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)

# Create a continuous norm to map from data points to colors
norm = plt.Normalize(dydx.min(), dydx.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
# Set the values used for colormapping
lc.set_array(dydx)
lc.set_linewidth(2)
line = axs[0].add_collection(lc)

plt.show()