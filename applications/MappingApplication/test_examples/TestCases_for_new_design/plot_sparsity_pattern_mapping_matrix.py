def spy(m, ax):
    from scipy.sparse import coo_matrix
    from matplotlib.patches import Rectangle
    if not isinstance(m, coo_matrix):
        m = coo_matrix(m)
    for (x, y) in zip(m.col, m.row):
        ax.add_artist(Rectangle(
            xy=(x-0.5, y-0.5), width=1, height=1))
    ax.set_xlim(-0.5, m.shape[1]-0.5)
    ax.set_ylim(-0.5, m.shape[0]-0.5)
    ax.invert_yaxis()
    ax.set_aspect(float(m.shape[0])/float(m.shape[1]))

import matplotlib.pylab as plt
import scipy.sparse
import scipy.io as sio
M = sio.mmread('MappingMatrixSerial')
# print(H)
fig, ax = plt.subplots()
plt.spy(M,precision=0.01, markersize=4)
plt.show()

a_x = sio.mmread('UpdateSystemVector')
# print(H)
fig, ax = plt.subplots()
plt.spy(a_x,precision=0.01, markersize=4)
plt.show()

b_x = sio.mmread('Update')
# print(H)
fig, ax = plt.subplots()
plt.spy(b_x,precision=0.01, markersize=4)
plt.show()

# spy(H,ax)   

# import matplotlib.pylab as pl
# import scipy.sparse as sps
# import scipy.io
# import sys
# A=scipy.io.mmread(sys.argv[1])
# pl.spy(A,precision=0.01, markersize=1)
# pl.show()