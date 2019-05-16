#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS17-18
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek

Author: philipp.bucher@tum.de, mate.pentek@tum.de

Description: Script for plotting aerodynamic forces and moments

Created on:  05.12.2015
Last update: 10.12.2017
'''
#===============================================================================
import sys
from matplotlib.pylab import *
import numpy as np

file2read = "FluidModelPart.Drag_structure_drag.dat" # file name of the drag results

# ======================================================

# *** read in the aerodynamic force/moment results ****************
simulTime = loadtxt(file2read, skiprows=3, usecols = (0,))
forceX = loadtxt(file2read, skiprows=3, usecols = (1,))
forceY = loadtxt(file2read, skiprows=3, usecols = (2,))
forceZ = loadtxt(file2read, skiprows=3, usecols = (3,))

try: # check if results for moments exist and read them if they exist
    momentX = loadtxt(file2read, skiprows=3, usecols = (4,))
    momentY = loadtxt(file2read, skiprows=3, usecols = (5,))
    momentZ = loadtxt(file2read, skiprows=3, usecols = (6,))
    momentsRead = True
    numberOfSubPlots = 2
except:
    momentsRead = False
    numberOfSubPlots = 1

# *** plot the aerodynamic forces/moments *************************
lw = 3 # line width in plot

plt.subplot(numberOfSubPlots,1,1) # plot the forces
plt.plot(simulTime, forceX, 'b-', linewidth=lw, label='Aerodyn. Force X without Ramp')
plt.plot(simulTime, forceY, 'r-', linewidth=lw, label='Aerodyn. Force Y without Ramp')
plt.plot(simulTime, forceZ, 'g-', linewidth=lw, label='Aerodyn. Force Z without Ramp')
plt.title('Aerodynamic Forces')
plt.ylabel('Force [N]')
plt.ticklabel_format(axis='y',style='sci',scilimits=(3,3)) # change y-axis numbers to scientific format
plt.grid(True)
plt.xlim(xmax=simulTime[-1]) # apply the maximum computation time as upper x-axis limit
plt.legend()
if not momentsRead: # plot time label if second plot doesn't exist
    plt.xlabel('Time [s]')

if momentsRead: # plot the moments if they have been read
    plt.subplot(numberOfSubPlots,1,2)
    plt.plot(simulTime, momentX, 'b-', linewidth=lw, label='Aerodyn. Moment X')
    plt.plot(simulTime, momentY, 'r-', linewidth=lw, label='Aerodyn. Moment Y')
    plt.plot(simulTime, momentZ, 'k-', linewidth=lw, label='Aerodyn. Moment Z')
    plt.title('Aerodynamic Moments')
    plt.xlabel('Time [s]')
    plt.ylabel('Moment [Nm]')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(3,3)) # change y-axis numbers to scientific format
    plt.grid(True)
    plt.xlim(xmax=simulTime[-1]) # apply the maximum computation time as upper x-axis limit
    plt.legend()

subplots_adjust(hspace=.35) # increase the space between the subplots
plt.show()