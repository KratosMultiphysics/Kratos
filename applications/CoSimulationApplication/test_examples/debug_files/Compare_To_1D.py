#Import of utilities
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.pyplot import xticks, yticks
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import sys
import os  # to be able to run Linux terminal commands
from matplotlib import cm
import warnings
import pickle
import matplotlib.animation as animation

# Setting of plot parameters
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

# Compare with or without inertia
inertia = True

# Note: the cases must have the same number of nodes

# Plot difference
diff = []
for i in range(1,101):
    if not inertia:
        file = f"../04_tube/CSM/Area_TS{i}"
    else:
        file = f"../04_tube_inert/CSM/Area_TS{i}"
    A = np.array(pd.read_csv(file, sep="\s+"))
    file = f"../04_tube2D/CSM/Area_TS{i}"
    B = np.array(pd.read_csv(file, sep="\s+"))
    diff.append(np.max(np.abs(B[:, 1] - A[:, 1])))

plt.figure()
plt.plot(A[:, 0], diff)
plt.title("Maximum value of absolute difference ifo z-coordinate")
plt.ylabel("differnce in area (m²)")
plt.xlabel("z-coordinate (m)")

plt.show()