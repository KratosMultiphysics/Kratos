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

matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

offset = 0.025

fig, ax = plt.subplots()
file = f"CFD2/Area_TS0"
A = np.array(pd.read_csv(file, sep="\s+", header=None))
x = A[:, 0] + offset
line, = ax.plot(x, A[:, 1], label='abaqus')
a0 = A[0, 1]
amp = 1e-7
plt.ylim([a0 - amp, a0 + amp])


file_inert = f"../04_tube/CSD_tangoparam/Area_TS0"
A_inert = np.array(pd.read_csv(file_inert, sep="\s+", header=None))
x_inert = A_inert[:, 0]
line_inert, = ax.plot(x_inert, A_inert[:, 1], label='1D with intertia')


def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(x))
    return line,

def init_inert():  # only required for blitting to give a clean slate.
    line_inert.set_ydata([np.nan] * len(x))
    return line_inert,

def animate(i):
    # print(i)
    file = f"CFD2/Area_TS{i}"
    A = np.array(pd.read_csv(file,sep="\s+", header=None))
    line.set_ydata(A[:,1])  # update the data.
    return line,

def animate_inert(i):
    # print(i)
    file_inert = f"../04_tube/CSD_tangoparam/Area_TS{i}"
    A_inert = np.array(pd.read_csv(file_inert,sep="\s+", header=None))
    line_inert.set_ydata(A_inert[:,1])  # update the data.
    return line_inert,

T = 100
ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=100, blit=False, save_count=50, repeat=True, frames=T)
ani_inert = animation.FuncAnimation(
    fig, animate_inert, init_func=init_inert, interval=100, blit=False, save_count=50, repeat=True, frames=T)

plt.plot(x, a0 * np.ones_like(A[:, 1]), color='k')
plt.legend()
plt.show()
