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

# Fluent Pipestructure
dir_pipstr = f"../fluent_pipeStructure/CSD/"
file = dir_pipstr + "Area_TS0"
A = np.array(pd.read_csv(file, sep="\s+", header=None))
x = A[:, 0]
line, = ax.plot(x, A[:, 1], label="2D Fluent + PipeStructure (CoCoNuT)")

a0 = A[0, 1]
amp = 5e-06
plt.ylim([a0 - amp / 4, a0 + amp])


def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(x))
    return line,


def animate(i):
    # print(i)
    file = dir_pipstr + f"Area_TS{i}"
    A = np.array(pd.read_csv(file, sep="\s+", header=None))
    line.set_ydata(A[:, 1])  # update the data.
    return line,


# 2D Fluent Abaqus CoCoNuT
dir_coco = '../tube2D_fluent_abaqus/CSM/'
file = dir_coco + f'CSM_Time0Surface0Nodes.dat'
A_coco = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
sort_id_coco = np.argsort(A_coco[:, 1])
pos0_coco = A_coco[sort_id_coco, :]
x_coco = pos0_coco[:, 1]
line_coco, = ax.plot(x_coco, A_coco[:, 0] ** 2 * np.pi, label="2D Fluent + Abaqus (CoCoNuT)")


def init_coco():  # only required for blitting to give a clean slate.
    line_coco.set_ydata([np.nan] * len(x_coco))
    return line_coco,


def animate_coco(i):
    # print(i)
    file = dir_coco + f'CSM_Time{i}Surface0Output.dat'
    A_coco = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))[sort_id_coco, :] + pos0_coco
    line_coco.set_xdata(A_coco[:, 1])
    line_coco.set_ydata(A_coco[:, 0] ** 2 * np.pi)  # update the data.
    return line_coco


# 2D Fluent Abaqus Tango
dir_Tango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube2D/FluentSolver0/'
file = dir_Tango + f'FSI1Refine0Time0Surface0NodesOrdered.dat'
A_tango = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
x_tango = A_tango[:, 1]
line_tango, = ax.plot(x_tango, A_tango[:, 2] ** 2 * np.pi, label="2D Fluent + Abaqus (Tango)", linestyle=':', color=line_coco.get_color())


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


# 3D Fluent Abaqus CoCoNuT
dir_3dcoco = '../tube3D_fluent_abaqus/CSM/'
file = dir_3dcoco + f'CSM_Time0Surface0Nodes.dat'
A_3dcoco = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
mask_3dcoco = (abs(A_3dcoco[:, 2]) < 1e-16) & (A_3dcoco[:, 1] > 0)
A_3dcoco = A_3dcoco[mask_3dcoco]
sort_id_3dcoco = np.argsort(A_3dcoco[:, 0])
pos0_3dcoco = A_3dcoco[sort_id_3dcoco, :]
print(pos0_3dcoco.shape)
x_3dcoco = pos0_3dcoco[:, 1]
line_3dcoco, = ax.plot(x_3dcoco, A_3dcoco[:, 1] ** 2 * np.pi, label="3D Fluent + Abaqus (CoCoNuT)")


def init_3dcoco():  # only required for blitting to give a clean slate.
    line_3dcoco.set_ydata([np.nan] * len(x_3dcoco))
    return line_3dcoco,


def animate_3dcoco(i):
    # print(i)
    file = dir_3dcoco + f'CSM_Time{i}Surface0Output.dat'
    A_3dcoco = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
    A_3dcoco = A_3dcoco[mask_3dcoco]
    A_3dcoco = A_3dcoco[sort_id_3dcoco, :] + pos0_3dcoco
    line_3dcoco.set_xdata(A_3dcoco[:, 0])
    line_3dcoco.set_ydata(A_3dcoco[:, 1] ** 2 * np.pi)  # update the data.
    return line_3dcoco,


# 3D Fluent Abaqus Tango
dir_3dtango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube3D/AbaqusSolver0/'
file = dir_3dtango + f'FSI1Refine0Time0Surface0Nodes.dat'
A_3dtango = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
mask_3dtango = (abs(A_3dtango[:, 2]) < 1e-16) & (A_3dtango[:, 1] > 0)
A_3dtango = A_3dtango[mask_3dtango]
sort_id_3dtango = np.argsort(A_3dtango[:, 0])
pos0_3dtango = A_3dtango[sort_id_3dtango, :]
print(pos0_3dtango.shape)
x_3dtango = pos0_3dtango[:, 0]
line_3dtango, = ax.plot(x_3dtango, A_3dtango[:, 2] ** 2 * np.pi, label="3D Fluent + Abaqus (Tango)", linestyle=':', color=line_3dcoco.get_color())


def init_3dtango():  # only required for blitting to give a clean slate.
    line_3dtango.set_ydata([np.nan] * len(x_3dtango))
    return line_3dtango,


def animate_3dtango(i):
    file = dir_3dtango + f'FSI1Refine0Time{i}Surface0Output.dat'
    A_3dtango = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
    A_3dtango = A_3dtango[mask_3dtango]
    A_3dtango = A_3dtango[sort_id_3dtango, :] + pos0_3dtango
    line_3dtango.set_xdata(A_3dtango[:, 0])
    line_3dtango.set_ydata(A_3dtango[:, 1] ** 2 * np.pi)  # update the data.
    return line_3dtango,


# Animation
T = 100
int = 100
ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=int, blit=False, save_count=50, repeat=True, frames=range(1, T))
ani_tango = animation.FuncAnimation(
    fig, animate_tango, init_func=init_tango, interval=int, blit=False, save_count=50, repeat=True, frames=range(1, T))
ani_coco = animation.FuncAnimation(
    fig, animate_coco, init_func=init_coco, interval=int, blit=False, save_count=50, repeat=True, frames=range(1, T))
ani_3dcoco = animation.FuncAnimation(
    fig, animate_3dcoco, init_func=init_3dcoco, interval=int, blit=False, save_count=50, repeat=True, frames=range(1, T))
ani_3dtango = animation.FuncAnimation(
    fig, animate_3dtango, init_func=init_3dtango, interval=int, blit=False, save_count=50, repeat=True, frames=range(1, T))

plt.plot(x, a0 * np.ones_like(A[:, 1]), color='k')
plt.legend()
plt.show()

