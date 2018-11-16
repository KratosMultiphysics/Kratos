# This code is designed to reproduce Fig. S7 from 2014 Vigolo

import numpy as np
import h5py
import matplotlib.pyplot as plt

L = 0.0048
z_plane_of_symmetry = -0.0024
n_intervals = 30
max_creation_time = 1.0
fluid_run_name = 'mesh_38288_nodes'
fluid_run_name = 'mesh_66576_nodes'
fluid_run_name = 'mesh_126189_nodes'
# fluid_run_name = 'mesh_242242_nodes'
snapshot_name = '1'
snapshot_name = '2'
a = 0.00024 / L
# snapshot_name = 't=0.04_in_box'

def IntervalIndex(x, bottom, top, number_of_intervals):
    if x < bottom:
        x = bottom
    elif x > top:
        x = top
        # raise ValueError('The coordinate z = ' + str(x) + ' is not within range!')
    return int(x / (top - bottom) * number_of_intervals)

def RobustQuotient(a, b):
    if a == b:
        return 1
    else:
        return a / b

def TooNew(creation_time):
    return creation_time > max_creation_time

def z_init(z):
    return abs(z - z_plane_of_symmetry) / L

def CalculateProbabilityOfGettingTrapped():
    half_width = abs(z_plane_of_symmetry)
    z_init_intervals = [[0.5 * i / n_intervals, 0.5 * (i + 1) / n_intervals] for i in range(n_intervals)]
    z_midpoints = [0.5 * (I[1] + I[0]) for I in z_init_intervals]
    histogram = np.zeros(len(z_init_intervals))
    trapped_hist = np.zeros(len(z_init_intervals))

    with h5py.File('all_particles.hdf5', 'r') as f:
        Ids = f['/Id'][:]
        Z0s = f['/Z0'][:]
        creation_times = f['/TIME'][:]
        too_new_to_include = [TooNew(time) for time in creation_times]
        ids_too_new_to_include = dict(zip(Ids, too_new_to_include))

        for z, Id in zip(Z0s, Ids):
            if not ids_too_new_to_include[Id]:
                histogram[IntervalIndex(z_init(z), 0., 0.5, n_intervals)] += 1

    with h5py.File('particles_snapshots.hdf5', 'r') as f:
        Ids_trapped = f['/' + fluid_run_name + '.hdf5/' + snapshot_name + '/Id'][:]
        Z0s_trapped = f['/' + fluid_run_name + '.hdf5/' + snapshot_name + '/Z0'][:]

        min_z =   float('inf')
        max_z = - float('inf')
        for Id, z in zip(Ids_trapped, Z0s_trapped):
            if not ids_too_new_to_include[Id]:
                z = z_init(z)
                min_z = min(min_z, z)
                max_z = max(max_z, z)
                trapped_hist[IntervalIndex(z , 0., 0.5, n_intervals)] += 1

    histogram = [RobustQuotient(h, H) for h, H in zip(trapped_hist, histogram)]
    plt.plot(z_midpoints, histogram, '.-')
    axes = plt.gca()
    axes.set_xlim([min_z, max_z])
    axes.set_ylim([0,1])
    plt.show()

CalculateProbabilityOfGettingTrapped()
