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

matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300

it = 1
TS = 1
# CFDfolder = "CFD-rb"
# CSMfolder = "CSM-rb"
CFDfolder = "CFD/"
CSMfolder = "CSM/"
Tangofolder = "/cfdfile2/data/fm/nicolas/Tango/cases/Data_Tube2D/"

### IMPORT DATA ###

# CoCoNuT Fluent
file_pos_CF = CFDfolder + f"nodes_update_timestep{TS}_thread3_Iter{it}.dat"  # xpos ypos id
file_pos2_CF = CFDfolder + f"nodes_update_timestep{TS}_thread3_Iter{it + 1}.dat"  # xpos ypos id
file_load_CF = CFDfolder + f"pressure_traction_timestep{TS}_thread3_Iter{it}.dat"  # WSSx WSSy p

file_nodes_CF = CFDfolder + f"nodes_timestep0_thread3.dat"  # x y id
file_faces_CF = CFDfolder + f"faces_timestep0_thread3.dat"  # x y id1 id2 (id3 id4)

CF_pos0 = np.array(pd.read_csv(file_nodes_CF, sep="\s+", header=None, skiprows=1))
CF_pos = np.array(pd.read_csv(file_pos_CF, sep="\s+", header=None, skiprows=1))
CF_pos2 = np.array(pd.read_csv(file_pos2_CF, sep="\s+", header=None, skiprows=1))
CF_lp = np.array(pd.read_csv(file_faces_CF, sep="\s+", header=None, skiprows=1))
CF_load = np.array(pd.read_csv(file_load_CF, sep="\s+", header=None, skiprows=1))

# CoCoNuT Abaqus
file_disp_CA = CSMfolder + f"CSM_Time{TS}Surface0Output_Iter{it}.dat"  # x y
file_load_CA = CSMfolder + f"CSM_Time{TS}Surface0Cpu0Input_Iter{it}.dat"  # p shear shear

file_nodes_CA = CSMfolder + "CSM_Time0Surface0Nodes.dat"
file_faces_CA = CSMfolder + "CSM_Time0Surface0Cpu0Faces.dat"  # elem_num lp_id y x

CA_disp = np.array(pd.read_csv(file_disp_CA, sep="\s+", header=None, skiprows=1))
CA_pos0 = np.array(pd.read_csv(file_nodes_CA, sep="\s+", header=None, skiprows=1))
CA_pos = CA_pos0 + CA_disp
file = CSMfolder + f"CSM_Time0Surface0Elements.dat"  # Tempory file for knowing how much to import from lp file
temp = np.array(pd.read_csv(file, header=None))
CA_lp = np.array(pd.read_csv(file_faces_CA, sep="\s+", header=None, skiprows=0))[:int(temp[0] * temp[1]), [2, 3]]
CA_load = np.array(pd.read_csv(file_load_CA, sep="\s+", header=None, skiprows=1))

# Tango Fluent
file_pos_TF = Tangofolder + f"FluentSolver0/FSI1Refine0Time{TS}Surface0Input_Iter{it}.dat"
file_pos2_TF = Tangofolder + f"FluentSolver0/FSI1Refine0Time{TS}Surface0Input_Iter{it + 1}.dat"
file_load_TF = Tangofolder + f"FluentSolver0/FSI1Refine0Time{TS}Surface0Output_Iter{it}.dat"  # x y id

file_nodes_TF = Tangofolder + f""  # x y
file_faces_TF = Tangofolder + f"FluentSolver0/FSI1Refine0Time0Surface0Faces.dat"  # x y id1 id2 (id3 id4)

TF_pos = np.array(pd.read_csv(file_pos_TF, header=None, skiprows=1, sep="\s+"))
TF_pos2 = np.array(pd.read_csv(file_pos2_TF, header=None, skiprows=1, sep="\s+"))
TF_lp = np.array(pd.read_csv(file_faces_TF, sep="\s+", header=None, skiprows=1))
TF_load = np.array(pd.read_csv(file_load_TF, sep="\s+", header=None, skiprows=1))

# Tango Abaqus
file_disp_TA = Tangofolder + f"AbaqusSolver0/FSI1Refine0Time{TS}Surface0Output_Iter{it}.dat"
file_load_TA = Tangofolder + f"AbaqusSolver0/FSI1Refine0Time{TS}Surface0Cpu0Input_Iter{it}.dat"

file_nodes_TA = Tangofolder + f"AbaqusSolver0/FSI1Refine0Time0Surface0Nodes.dat"
file_faces_TA = Tangofolder + f"AbaqusSolver0/FSI1Refine0Time0Surface0Faces.dat"  # elem_num lp_id x y

TA_disp = np.array(pd.read_csv(file_disp_TA, header=None, skiprows=1, sep="\s+"))
TA_pos0 = np.array(pd.read_csv(file_nodes_TA, header=None, skiprows=1, sep="\s+"))
TA_pos = TA_pos0 + TA_disp
TA_lp = np.array(pd.read_csv(file_faces_TA, sep="\s+", header=None, skiprows=1))[:, [2, 3]]
TA_load = np.array(pd.read_csv(file_load_TA, sep="\s+", header=None, skiprows=1))


### SORT ###


def sort(ar, index):
    if type(index) == int:
        sort_id = np.argsort(ar[:, index])
    else:
        sort_id = np.argsort(index)
    tmp = ar[sort_id, :]
    return tmp


# lsF = [CF_pos0, CF_pos, CF_pos2, TF_pos, TF_pos2, TA_pos0, TA_pos, TA_disp]
# lsA = [CA_pos0, CA_pos, CA_disp]
CF_pos0 = sort(CF_pos0, 0)
CF_pos = sort(CF_pos, 0)
CF_pos2 = sort(CF_pos2, 0)
TF_pos = sort(TF_pos, 0)
TF_pos2 = sort(TF_pos2, 0)
TA_pos0 = sort(TA_pos0, 0)
TA_pos = sort(TA_pos, 0)
TA_disp = sort(TA_disp, 0)
CA_pos0 = sort(CA_pos0, 1)
CA_pos = sort(CA_pos, 1)
CA_disp = sort(CA_disp, 1)

# lsload = [CF_load, TF_load, CF_load]
# lslp = [CF_lp, TF_lp, CF_lp]
CF_load = sort(CF_load, CF_lp[:, 0])
TF_load = sort(TF_load, TF_lp[:, 0])
TA_load = sort(TA_load, TA_lp[:, 0])
CA_load = sort(CA_load, CA_lp[:, 1])
CF_lp = sort(CF_lp, 0)
TF_lp = sort(TF_lp, 0)
TA_lp = sort(TA_lp, 0)
CA_lp = sort(CA_lp, 1)


### POSITION ###

plt.figure(1)

# COCO Fluent position
plt.plot(CF_pos0[:, 0], CF_pos0[:, 1], label="Initial")
plt.plot(CF_pos[:, 0], CF_pos[:, 1], label=f"CoCo_x{it}")
plt.plot(CF_pos2[:, 0], CF_pos2[:, 1], label=f"CoCo_x{it + 1}")

# COCO Abaqus position
plt.plot(CA_pos[:, 1], CA_pos[:, 0], label=f"Coco_xt{it}")

# TANGO Fluent position

plt.plot(TF_pos[:, 0], TF_pos[:, 1], ':', label=f"Tango_x{it}")
plt.plot(TF_pos2[:, 0], TF_pos2[:, 1], ':', label=f"Tango_x{it + 1}")

# TANGO Abaqus position

plt.plot(TA_pos[:, 0], TA_pos[:, 1], ':', label=f"Tango_xt{it}")

plt.title('Position')
plt.legend()

### LOAD ###

# COCO load
plt.figure(2)
plt.plot(CA_lp[:, 1], CA_load[:, 0], '--', label=f"CoCo_p{it}", marker='o')
plt.plot(CF_lp[:, 0], CF_load[:, 2], label=f"CoCo_pt{it}", marker='o')
plt.ylabel("pressure")
plt.title("pressure")

plt.figure(3)
plt.plot(CA_lp[:, 1], CA_load[:, 2], '--', label=f"CoCo_p{it}")
plt.plot(CF_lp[:, 0], CF_load[:, 0], label=f"CoCo_pt{it}")
plt.ylabel("axial shear")
plt.title("axial shear")

plt.figure(4)
plt.plot(CA_lp[:, 1], CA_load[:, 1], '--', label=f"CoCo_p{it}")
plt.plot(CF_lp[:, 0], CF_load[:, 1], label=f"CoCo_pt{it}")
plt.ylabel("radial shear")
plt.title("radial shear")

# TANGO load
if True:
    plt.figure(2)
    plt.plot(TA_lp[:, 0], TA_load[:, 0], ':', label=f"Tango_p{it}", marker='o')
    plt.plot(TF_lp[:, 0], TF_load[:, 0], ':', label=f"Tango_pt{it}", marker='o')
    plt.legend()

    plt.figure(3)
    plt.plot(TA_lp[:, 0], TA_load[:, 1], ':', label=f"Tango_p{it}")
    plt.plot(TF_lp[:, 0], TF_load[:, 1], ':', label=f"Tango_pt{it}")
    plt.legend()

    plt.figure(4)
    plt.plot(TA_lp[:, 0], TA_load[:, 2], ':', label=f"Tango_p{it}")
    plt.plot(TF_lp[:, 0], TF_load[:, 2], ':', label=f"Tango_pt{it}")
    plt.legend()

    plt.figure(5)
    plt.plot(TA_lp[:, 0], TA_load[:, 0] - CA_load[:, 0], label="diff_load_Abaqus")

    # print(TF_load[:,0].shape, TA_load[:,0].shape)

    CTx = TF_load[:, 0]
    CTy = TF_load[:, 1]
    CTz = TF_load[:, 2]

    CCx = CF_load[:, 2]
    CCy = CF_load[:, 0]
    CCz = CF_load[:, 1]

    plt.plot(TF_lp[:, 0], CTx - CCx, label="diff_pressure_Fluent")
    plt.plot(TF_lp[:, 0], CTy - CCy, label="diff_WSSx_Fluent")
    plt.plot(TF_lp[:, 0], CTz - CCz, label="diff_WSSy_Fluent")

    plt.title('difference')
    plt.legend()

plt.show()
