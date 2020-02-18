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

fig, ax = plt.subplots()
dir_pipstr = f"../fluent_pipeStructure/CSD/"
file = dir_pipstr + "Area_TS0"
A = np.array(pd.read_csv(file, sep="\s+", header=None))
x = A[:, 0]
line, = ax.plot(x, A[:, 1], label="fluent_pipestructure")
a0 = A[0, 1]
amp = 5e-06
plt.ylim([a0 - amp, a0 + amp])

def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(x))
    return line,


def animate(i):
    # print(i)
    file = dir_pipstr + f"Area_TS{i}"
    A = np.array(pd.read_csv(file, sep="\s+", header=None))
    line.set_ydata(A[:, 1])  # update the data.
    return line,


dir_Tango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube2D/FluentSolver0/'
file = dir_Tango + f'FSI1Refine0Time0Surface0NodesOrdered.dat'
A_tango = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
x_tango = A_tango[:, 1]
line_tango, = ax.plot(x_tango, A_tango[:, 2] ** 2 * np.pi, label="Tango")


def init_tango():  # only required for blitting to give a clean slate.
    line_tango.set_ydata([np.nan] * len(x_tango))
    return line_tango,


def animate_tango(i):
    # print(i)
    file = dir_Tango + f'FSI1Refine0Time{i}Surface0NodesOrdered.dat'
    A_tango = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
    line_coco.set_xdata(A_tango[:, 0])
    line_tango.set_ydata(A_tango[:, 1] ** 2 * np.pi)  # update the data.
    return line_tango,

dir_coco = 'CSM/'
file = dir_coco + f'CSM_Time0Surface0Nodes.dat'
A_coco = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
sort_id = np.argsort(A_coco[:, 1])
pos0_coco = A_coco[sort_id, :]
x_coco = pos0_coco[:, 1]
line_coco, = ax.plot(x_coco, A_coco[:, 0] ** 2 * np.pi, label="coco")


def init_coco():  # only required for blitting to give a clean slate.
    line_coco.set_ydata([np.nan] * len(x_coco))
    return line_coco,


def animate_coco(i):
    # print(i)
    file = dir_coco + f'CSM_Time{i}Surface0Output.dat'
    A_coco = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))[sort_id, :] + pos0_coco
    line_coco.set_xdata(A_coco[:, 1])
    line_coco.set_ydata(A_coco[:, 0] ** 2 * np.pi)  # update the data.
    return line_coco,

T = 100
ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=100, blit=False, save_count=50, repeat=True, frames=range(1, T))
ani_tango = animation.FuncAnimation(
    fig, animate_tango, init_func=init_tango, interval=100, blit=False, save_count=50, repeat=True, frames=range(1, T))
ani_coco = animation.FuncAnimation(
    fig, animate_coco, init_func=init_coco, interval=100, blit=False, save_count=50, repeat=True, frames=range(1, T))

plt.plot(x, a0 * np.ones_like(A[:, 1]), color='k')
plt.legend()
plt.show()

