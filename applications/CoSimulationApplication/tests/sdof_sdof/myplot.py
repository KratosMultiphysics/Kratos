import sys
from matplotlib.pylab import *
import numpy as np

file2read = "fsi_sdof/fsi_sdof_cfd_results_force.dat" # file name of the drag results

# ======================================================

# *** read in the aerodynamic force/moment results ****************
simulTime = loadtxt(file2read, skiprows=3, usecols = (0,))
forceX = loadtxt(file2read, skiprows=3, usecols = (1,))
forceY = loadtxt(file2read, skiprows=3, usecols = (2,))
forceZ = loadtxt(file2read, skiprows=3, usecols = (3,))

# *** plot the aerodynamic forces/moments *************************
lw = 1

plt.figure(1) # plot the forces
plt.plot(simulTime, forceX, 'b-', linewidth=lw, label='Aerodyn. Force X without Ramp')
plt.plot(simulTime, forceY, 'r-', linewidth=lw, label='Aerodyn. Force Y without Ramp')
plt.plot(simulTime, forceZ, 'g-', linewidth=lw, label='Aerodyn. Force Z without Ramp')
plt.title('Aerodynamic Forces')
plt.ylabel('Force [N]')
plt.ticklabel_format(axis='y',style='sci',scilimits=(3,3)) # change y-axis numbers to scientific format
plt.grid(True)
plt.xlim(xmax=simulTime[-1]) # apply the maximum computation time as upper x-axis limit
plt.legend()

subplots_adjust(hspace=.35) # increase the space between the subplots
plt.savefig('fsi_sdof/myfig_cfd.png')
#plt.show()

file2read = "fsi_sdof/fsi_sdof_cfd_results_disp.dat" # file name of the drag results

# ======================================================

# *** read in the aerodynamic force/moment results ****************
simulTime = loadtxt(file2read, skiprows=3, usecols = (0,))
forceX = loadtxt(file2read, skiprows=3, usecols = (1,))
# forceY = loadtxt(file2read, skiprows=3, usecols = (2,))
# forceZ = loadtxt(file2read, skiprows=3, usecols = (3,))

# *** plot the aerodynamic forces/moments *************************
lw = 1

plt.figure(2) # plot the forces
plt.plot(simulTime, forceX, 'b-', linewidth=lw, label='Aerodyn. Force X without Ramp')
# plt.plot(simulTime, forceY, 'r-', linewidth=lw, label='Aerodyn. Force Y without Ramp')
# plt.plot(simulTime, forceZ, 'g-', linewidth=lw, label='Aerodyn. Force Z without Ramp')
plt.title('Aerodynamic Forces')
plt.ylabel('Force [N]')
plt.ticklabel_format(axis='y',style='sci',scilimits=(3,3)) # change y-axis numbers to scientific format
plt.grid(True)
plt.xlim(xmax=simulTime[-1]) # apply the maximum computation time as upper x-axis limit
plt.legend()


subplots_adjust(hspace=.35) # increase the space between the subplots
plt.savefig('fsi_sdof/myfig_sdof.png')
#plt.show()