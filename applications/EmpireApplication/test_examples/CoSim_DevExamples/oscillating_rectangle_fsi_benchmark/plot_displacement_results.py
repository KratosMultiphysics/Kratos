#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS17-18
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek

Author: philipp.bucher@tum.de, mate.pentek@tum.de

Description: Script for plotting data over time

Created on:  05.12.2015
Last update: 10.12.2017
'''
#===============================================================================

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from pylab import *

# ======================================================
import json

file2read = "results_sdof.dat"
# ======================================================

# read the file
simulTime = loadtxt(file2read, skiprows=1, usecols = (0,))
timeData = loadtxt(file2read, skiprows=1, usecols = (1,))

# set up the plot
fig = plt.figure()
ax = plt.axes(xlim=(min(simulTime), max(simulTime)), ylim=(min(timeData), max(timeData)))
line, = ax.plot([], [], lw=2)
ax.set_xlabel("Time [s]")
ax.set_ylabel("Displacement [m]")
ax.set_title('Displacement Results for SDoF')
ax.plot(simulTime, timeData, "-b", lw=1.0)
plt.ticklabel_format(axis='y',style='sci',scilimits=(3,3)) # change y-axis numbers to scientific format
plt.grid(True)

plt.show()
