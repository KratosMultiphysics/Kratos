import matplotlib.pyplot as plt
from numpy import loadtxt

file2read = 'results_sdof_ref.dat'

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