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

directory = './CFD-rb/'
dir_abaq = './CSM-rb/'

file = dir_abaq+"CSM_Time0Surface0Nodes.dat"
CSM_pos0 = np.array(pd.read_csv(file,header=None,skiprows=1,sep="\s+"))

dir_Tango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube2D/FluentSolver0/'
T = [i for i in range(1,2)]

plt.figure()
for t in T:
    file = directory+f'nodes_update_timestep{t}_thread3.dat'
    pos = np.array(pd.read_csv(file,header=None,skiprows=1,sep="\s+"))
    sort_id = np.argsort(pos[:,0])
    pos=pos[sort_id,:]
    p = plt.plot(pos[:,0],pos[:,1],marker='o',label=f"CoCoNuT_{t}")

    file = dir_Tango+f'FSI1Refine0Time{t}Surface0NodesOrdered.dat'
    posT = np.array(pd.read_csv(file,header=None,skiprows=1,sep="\s+"))
    plt.plot(posT[:,0],posT[:,1],linestyle='--', marker='x',label=f"Tango_{t}",color= p[0].get_color())

    file = dir_abaq+f"CSM_Time{t}Surface0Output.dat"
    posA = CSM_pos0 + np.array(pd.read_csv(file,header=None,skiprows=1,sep="\s+"))
    sort_id = np.argsort(posA[:,1])
    posA=posA[sort_id,:]
    plt.plot(posA[:, 1], posA[:, 0], linestyle=':', marker='x', label=f"Abaqus_{t}", color= p[0].get_color())

plt.legend(ncol=2)

plt.figure()
t = 1
iter = [1,2]
for i in iter:
    file = dir_abaq+f"CSM_Time{t}Surface0Output_Iter{i}.dat"
    posA = CSM_pos0 + np.array(pd.read_csv(file,header=None,skiprows=1,sep="\s+"))
    sort_id = np.argsort(posA[:,1])
    posA=posA[sort_id,:]
    plt.plot(posA[:, 1], posA[:, 0], linestyle=':', marker='x', label=f"Abaqus_Iter{i}", color= p[0].get_color())
    plt.legend()


plt.show()