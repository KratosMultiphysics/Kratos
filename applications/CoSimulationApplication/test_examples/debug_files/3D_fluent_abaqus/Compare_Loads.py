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



file_loadTango = "/cfdfile2/data/fm/lucas/Tango_Test/Data/AbaqusSolver0/FSI1Refine0Time1Surface0Cpu0Input.dat"
file_coordTango = "/cfdfile2/data/fm/lucas/Tango_Test/Data/AbaqusSolver0/FSI1Refine0Time0Surface0Faces.dat"

coordTango = np.array(pd.read_csv(file_coordTango, sep="\s+"))
loadTango = np.array(pd.read_csv(file_loadTango,  sep="\s+",header = None, skiprows=1))
n_Tango = loadTango.shape[0]

file_loadCoCo =  "/cfdfile2/data/fm/lucas/Kratos/applications/CoSimulationApplication/test_examples/tube3D_fluent_abaqus/CSM-rb/CSM_Time1Surface0Cpu0Input.dat"
file_coordCoCo = "/cfdfile2/data/fm/lucas/Kratos/applications/CoSimulationApplication/test_examples/tube3D_fluent_abaqus/CSM-rb/CSM_Time0Surface0CFaces.dat"

coordCoCo = np.array(pd.read_csv(file_coordTango,  sep="\s+"))
loadCoCo = np.array(pd.read_csv(file_loadCoCo,  sep="\s+", header = None, skiprows = 1))
n_coco = loadCoCo.shape[0]

plt.figure()
plt.scatter(coordTango[:,2],loadTango[:,0],label="pressure Tango")
plt.scatter(coordCoCo[:,2],loadCoCo[:,0],label="pressure CoCo")
plt.legend()

plt.figure()
plt.scatter(coordCoCo[:,2],loadCoCo[:,0]-loadTango[:,0],label="pressure CoCo-Tango")
plt.legend()

plt.show()