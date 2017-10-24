# This code is designed to reproduce Fig. S7 from 2014 Vigolo

import numpy as np
import h5py
import matplotlib.pyplot as plt

z_plane_of_symmetry = -0.0024
n_intervals = 10
fluid_run_name = 'mesh_38288_nodes'
snapshot_name = 't=0.2900000000000001_in_box'

def IntervalIndex(x, bottom, top, number_of_intervals):
    if x < bottom:
        x = bottom
    elif x > top:
        x = top
        # raise ValueError('The coordinate z = ' + str(x) + ' is not within range!')
    return int(x / (top - bottom) * number_of_intervals)

def CalculateProbabilityOfGettingTrapped():
    z_min = z_plane_of_symmetry
    z_init_intervals = [[z_min + i       * abs(2 * z_min) / n_intervals,
                         z_min + (i + 1) * abs(2 * z_min) / n_intervals] for i in range(n_intervals)]
    z_midpoints = [0.5 * (i[1] + i[0]) for i in z_init_intervals]

    with h5py.File('all_particles.hdf5', 'r') as f:
        Z0s = f['/Z0'][:]
        N = len(Z0s)
        histogram = [0.0 for i in z_init_intervals]
        for z in Z0s:
            histogram[IntervalIndex(z, z_plane_of_symmetry, 0., n_intervals)] += 1

    with h5py.File('particles_snapshot.hdf5', 'r') as f:
        Z0s_trapped = f['/' + fluid_run_name + '.hdf5/' + snapshot_name + '/Z0'][:]
        trapped_hist = [0.0 for i in z_init_intervals]
        for z in Z0s_trapped:
            trapped_hist[IntervalIndex(z, z_plane_of_symmetry, 0., n_intervals)] += 1
        N_trapped = len(Z0s_trapped)

    def RobustQuotient(a, b):
        if a == b:
            return 1
        else:
            return a / b
    histogram = [RobustQuotient(h, H) for h, H in zip(trapped_hist, histogram)]
    plt.plot(z_midpoints, histogram, '.')
    plt.show()

CalculateProbabilityOfGettingTrapped()
