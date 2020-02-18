# Import of utilities
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

# Test_examples to plot
plot_04_tube = True
plot_04_tube_inert = True
plot_04_tube2D = True

# Initialization
fig, ax = plt.subplots()
def file(i):
    return f"/CSM/Area_TS{i}"
path = "../"

# Plot definitions
if plot_04_tube:
    name = "04_tube"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube = A[:, 0]
    line_04_tube, = ax.plot(x_04_tube, A[:, 1], label=name)

    def init_04_tube():  # only required for blitting to give a clean slate.
        line_04_tube.set_ydata([np.nan] * len(x_04_tube))
        return line_04_tube,

    def animate_04_tube(i):
        tmp = path + "04_tube" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube.set_ydata(A[:, 1])  # update the data.
        return line_04_tube,

if plot_04_tube_inert:
    name = "04_tube_inert"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube_inert = A[:, 0]
    line_04_tube_inert, = ax.plot(x_04_tube_inert, A[:, 1], label=name)

    def init_04_tube_inert():  # only required for blitting to give a clean slate.
        line_04_tube_inert.set_ydata([np.nan] * len(x_04_tube_inert))
        return line_04_tube_inert,

    def animate_04_tube_inert(i):
        tmp = path + "04_tube_inert" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube_inert.set_ydata(A[:, 1])  # update the data.
        return line_04_tube_inert,

if plot_04_tube2D:
    name = "04_tube2D"
    tmp = path + name + file(0)
    A = np.array(pd.read_csv(tmp, sep="\s+"))
    x = x_04_tube2D = A[:, 0]
    line_04_tube2D, = ax.plot(x_04_tube2D, A[:, 1], label=name)

    def init_04_tube2D():  # only required for blitting to give a clean slate.
        line_04_tube2D.set_ydata([np.nan] * len(x_04_tube2D))
        return line_04_tube2D,

    def animate_04_tube2D(i):
        tmp = path + "04_tube2D" + file(i)
        A = np.array(pd.read_csv(tmp, sep="\s+"))
        line_04_tube2D.set_ydata(A[:, 1])  # update the data.
        return line_04_tube2D,

# Setting scale
a0 = A[0, 1]
amp = 3e-07
plt.ylim([a0 - amp, a0 + amp])

# Animation parameters
interv = 100  # Interval between frames in ms
T = 100  # Number of frames

# Animations
if plot_04_tube:
    ani_04_tube = animation.FuncAnimation(fig, animate_04_tube, init_func=init_04_tube, interval=interv, blit=False, save_count=50, repeat=True, frames=T)
if plot_04_tube_inert:
    ani_04_tube_inert = animation.FuncAnimation(fig, animate_04_tube_inert, init_func=init_04_tube_inert, interval=interv, blit=False, save_count=50, repeat=True, frames=T)
if plot_04_tube2D:
    ani_04_tube2D = animation.FuncAnimation(fig, animate_04_tube2D, init_func=init_04_tube2D, interval=interv, blit=False, save_count=50, repeat=True, frames=T)

plt.plot(x, A[:, 1], color='k')
plt.title("Animation of the area of each crosssection")
plt.ylabel("area (m²)")
plt.xlabel("z-coordinate (m)")
plt.legend()
plt.show()