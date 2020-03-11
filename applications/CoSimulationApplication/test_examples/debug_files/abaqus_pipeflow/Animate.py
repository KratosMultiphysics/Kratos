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

matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

fig, ax = plt.subplots()
file = f"CSM/CSM_Time0Surface0Nodes.dat"
A = np.array(pd.read_csv(file, sep="\s+",header=None,skiprows=1))
x = A[:,1]
y0 = A[:,0]*0.0
isort = np.argsort(x)
line, = ax.plot(x[isort], y0[isort],marker='x')
ax.set_ylim([0,1e-05])

T_end = 101


def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(x))
    return line,


def animate(i):
    # print(i)
    file = f"CSM/CSM_Time{i}Surface0Output.dat"
    A = np.array(pd.read_csv(file,sep="\s+",header=None,skiprows=1))
    y=y0+A[:,0]
    line.set_ydata(y[isort])  # update the data.
    return line,

fig2, ax2 = plt.subplots()
file2 = f"CSM/CSM_Time0Surface0Cpu0Faces.dat"
B = np.array(pd.read_csv(file2, sep="\s+",header=None))
x2 = B[:150,3]
l0 = B[:150,0]*0.0
isort2 = np.argsort(x2)
line2, = ax2.plot(x2[isort2], l0[isort2],marker='x')
ax2.set_ylim([0,100])

def init2():  # only required for blitting to give a clean slate.
    line2.set_ydata([np.nan] * len(x))
    return line2,


def animate2(i):
    # print(i)
    file2 = f"CSM/CSM_Time{i}Surface0Cpu0Input.dat"
    B = np.array(pd.read_csv(file2,sep="\s+",header=None,skiprows=1))
    y=B[:,0]
    line2.set_ydata(y[isort2])  # update the data.
    return line2,
#
anims=[]

anims.append(animation.FuncAnimation(fig, animate, init_func=init, interval=25, blit=True, save_count=50,repeat=True,frames=range(1,T_end)))
anims.append(animation.FuncAnimation(fig2, animate2, init_func=init, interval=25, blit=True, save_count=50,repeat=True,frames=range(1,T_end)))

### To get an mp4-file
#matplotlib.verbose.set_level("helpful")
# print(plt.rcParams['animation.ffmpeg_path'])
# plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'
# print(plt.rcParams['animation.ffmpeg_path'])
# writer = animation.FFMpegFileWriter(codec='mpeg1video', fps=24, bitrate=2000)

#
# anims[0].save("Displacement.mpeg", writer=writer)
# anims[1].save("Pressure.mpeg", writer=writer)



plt.show()