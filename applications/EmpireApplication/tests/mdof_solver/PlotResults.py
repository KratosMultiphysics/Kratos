import matplotlib.pyplot as plt
from numpy import loadtxt

#file2read = 'results_mdof_sdof_ref.dat'
#file2read = 'results_mdof_cantilever_shear_2d_ref.dat'
#file2read = 'results_mdof_bridge_2dof_ref.dat'
file2read = 'results_mdof_generic_ref.dat'

col_counter = 0
X = loadtxt(file2read, skiprows=1, usecols = (col_counter,))
col_counter += 1
try:
  while True:
    Y = loadtxt(file2read, skiprows=1, usecols = (col_counter,))
    col_counter += 1
except:
  print('reached last column')

plt.plot(X, Y)
plt.show()
