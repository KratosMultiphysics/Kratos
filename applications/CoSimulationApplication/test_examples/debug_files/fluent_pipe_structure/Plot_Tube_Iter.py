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

dir_Tango = '/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube2D/'
dir_CoCo = '../fluent_pipe_structure/'

T = [60]
it_max = 10
CFDfolder = "CFD/"
CSMfolder = "CSM/"
Tango = True

if Tango:
    file = dir_Tango+"AbaqusSolver0/"+f"FSI1Refine0Time0Surface0Nodes.dat"
    pos0_Tango = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
else:
    dir_Tango = ''

# file = dir_CoCo + CSMfolder + f"CSM_Time0Surface0Nodes.dat"
# pos0_CoCo = np.array(pd.read_csv(file, header=None, skiprows=1, sep="\s+"))
# print(os.getcwd())
# print()

for t in T:
    for i in range(1, it_max):
        file_TF = dir_Tango+"FluentSolver0/"+f"FSI1Refine0Time{t}Surface0Input_Iter{i}.dat"
        file_TA = dir_Tango+"AbaqusSolver0/"+f"FSI1Refine0Time{t}Surface0Output_Iter{i}.dat"

        file_CF = dir_CoCo + CFDfolder + f"nodes_update_timestep{t}_thread3_Iter{i}.dat"
        file_CA = dir_CoCo + CSMfolder + f"CSM_Time{t}Surface0Output_Iter{i}.dat"
        b = 0
        b2 = 0

        if os.path.exists(file_TF):
            pos = np.array(pd.read_csv(file_TF, header=None, skiprows=1, sep="\s+"))
            plt.figure(4 * t + 1)
            sort_id = np.argsort(pos[:, 0])
            pos = pos[sort_id, :]
            p = plt.plot(pos[:, 0], pos[:, 1], label=f"Tango_x{i}")
            b = 1

        if os.path.exists(file_CF):
            pos = np.array(pd.read_csv(file_CF, header=None, skiprows=1, sep="\s+"))
            plt.figure(4 * t + 1)
            sort_id = np.argsort(pos[:, 0])
            pos = pos[sort_id, :]
            if b == 1:
                plt.plot(pos[:, 0], pos[:, 1], linestyle='--', label=f"CoCoNuT_x{i}", color=p[0].get_color())
            else:
                plt.plot(pos[:, 0], pos[:, 1], linestyle='-', label=f"CoCoNuT_x{i}")

        if os.path.exists(file_TA):
            pos = np.array(pd.read_csv(file_TA, header=None, skiprows=1, sep="\s+")) + pos0_Tango
            plt.figure(4 * t + 2)
            # pos=pos[:,[1,0]]
            sort_id = np.argsort(pos[:, 0])
            pos = pos[sort_id,:]
            p = plt.plot(pos[:,0],pos[:,1],label=f"Tango_xt{i}")
            b2 = 1

        if os.path.exists(file_CA):
            pos = np.array(pd.read_csv(file_CA, header=None, skiprows=1, sep="\s+"))
            plt.figure(4 * t + 2)
            # pos=pos[:,[1,0]]
            sort_id = np.argsort(pos[:, 0])
            pos = pos[sort_id, :]
            if b2 == 1:
                plt.plot(pos[:, 0], pos[:, 1], linestyle=':', label=f"CoCoNuT_xt{i}", color=p[0].get_color())
            else:
                plt.plot(pos[:, 0], pos[:, 1], linestyle=':', label=f"CoCoNuT_xt{i}")
        plt.figure(4 * t + 1)
        plt.title(f"Timestep {t}")
        leg = plt.legend()
        # leg.draggable()
        plt.figure(4 * t + 2)
        leg = plt.legend()
        # leg.draggable()

plt.show()