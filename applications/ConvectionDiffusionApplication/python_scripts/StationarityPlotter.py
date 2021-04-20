# # Import system python
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def GetFieldDataFile(mdpa_name):
    model_part_name = mdpa_name
    hf = h5py.File(model_part_name, 'r')
    velocity_error_vector = []
    concentration_error_vector = []
    group_names = list(hf.keys())
    group_names = list(map(int,group_names))
    group_names.sort(reverse=False)
    times_str = list([str(key) for key in hf.keys() if 'time' in hf['/' + key].attrs])
    times = np.array([float(hf[key].attrs['time']) for key in times_str])
    times = np.sort(times)
    for i in (group_names):
        velocity_error_vector = np.append(velocity_error_vector, hf[str(group_names[i-1])+'/VELOCITY_ERROR'])
        concentration_error_vector = np.append(concentration_error_vector, hf[str(group_names[i-1])+'/CONCENTRATION_ERROR'])
        
    return velocity_error_vector, concentration_error_vector, times


def PlotMeshConvergence(velocity, concentration, time):
    
    plt.subplot(1, 2, 1)
    plt.plot(time, velocity, color='b', marker= 'o')
    plt.title("Velocity", fontsize=17)
    plt.grid()
    # plt.xlim(0, 0.4)
    plt.yscale('log')
    plt.ylabel('L2  error  norm', fontsize=15, labelpad=5)
    plt.xlabel('time [s]', fontsize=15, labelpad=5)
    plt.tick_params(axis='both', which='major', labelsize=15)
    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)
    # plt.suptitle("Stationarity", fontsize=20)
    # figure = plt.gcf() # get current figure
    # plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
    # figure.set_size_inches(14, 11)
    # plt.savefig('L2 error norm' + '.pdf', format='pdf', dpi=1000)
    # plt.show()

    # plt.subplot(1, 2, 2)
    # plt.plot(step, concentration, color='b', marker='o', label = '$Concentration$')
    # plt.title("Concentration", fontsize=17)
    # plt.grid()
    # plt.xlim(0, 0.75)
    # plt.yscale('log')
    # # plt.ylabel('L2  error  norm', fontsize=20, labelpad=25)
    # plt.xlabel('steps', fontsize=15, labelpad=5)
    # # lgnd = plt.legend(loc = 'lower left', prop={'size':15}, frameon=False)
    # plt.tick_params(axis='both', which='major', labelsize=15)
    # # plt.semilogy()
    # ax = plt.gca()
    # ax.tick_params(axis='x', pad=20)
    # ax.tick_params(axis='y', pad=10)
    
    # plt.subplot(1, 2, 2)
    # plt.plot(time, width, color='b', marker= 'o')
    # plt.title("Diffusion layer", fontsize=17)
    # plt.grid()
    # # plt.xlim(0, 0.4)
    # plt.yscale('log')
    # plt.ylabel('Width', fontsize=15, labelpad=5)
    # plt.xlabel('time [s]', fontsize=15, labelpad=5)
    # plt.tick_params(axis='both', which='major', labelsize=15)
    # ax = plt.gca()
    # ax.tick_params(axis='x', pad=20)
    # ax.tick_params(axis='y', pad=10)

    plt.suptitle("Stationarity", fontsize=20)
    figure = plt.gcf() # get current figure
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
    figure.set_size_inches(14, 11)
    plt.savefig('L2 error norm' + '.pdf', format='pdf', dpi=1000)
    plt.show()
        
d = {}
velocity_error_v = []
concentration_error_v = []
times_v = []
width_v = []
d['mesh_4'+'_dirPath'] = '/home/aitor/Escritorio/Norouzi_mesh_4.gid/Norouzi_mesh_4.hdf5'
velocity_error_v, concentration_error_v, times_v = GetFieldDataFile(d['mesh_4'+'_dirPath'])
PlotMeshConvergence(velocity_error_v, concentration_error_v, times_v)

            




        