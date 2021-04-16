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

def GetWidth(mdpa_name):
    model_part_name = mdpa_name
    hf = h5py.File(model_part_name, 'r')
    target_node = 5402
    coordinates = []
    nodes_id = []
    velocity = []
    search_id = []
    concentration = []
    found_concentration = []
    found_coordinate_x = []
    found_coordinates = []
    node_ids = []
    width = []
    group_names = list(hf.keys())
    group_names = list(map(int,group_names))
    group_names.sort(reverse=False)
    for i in (group_names):
        n_size = float(hf['/' + str(i)].attrs['element_size'])
        coordinates = hf[str(group_names[i-1])+'/COORDINATES']
        nodes_id = hf[str(group_names[i-1])+'/ID']w
        # velocity = hf[str(group_names[i-1])+'/VELOCITIES']
        concentration = hf[str(group_names[i-1])+'/TEMPERATURES']
        target_coordinate = coordinates[target_node-1]
        target_tolerance = [target_coordinate[0]-n_size, target_coordinate[0]+n_size]
        for node_i, node in enumerate(coordinates):
            test_value = coordinates[node_i][0]
            if  test_value > target_tolerance[0] and test_value < target_tolerance[1]:
                # found_coordinates = coordinates[node_i]
                # found_coordinate_x = np.append(found_coordinate_x, found_coordinates[0])
                found_concentration = np.append(found_concentration, concentration[node_i])
                node_ids = np.append(node_ids, node_i)
        c0 = np.amax(found_concentration)
        c1 = np.amin(found_concentration)
        concentration_delta = [1.05*c1, 0.05*c0]
        diameter = c0-c1
        delta_1 = np.abs(c0-concentration_delta[1])
        delta_2 = np.abs(c1-concentration_delta[0])
        width = np.append(width, diameter-(delta_1+delta_2))
    return width

def PlotMeshConvergence(velocity, concentration, time, width):
    
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
    
    plt.subplot(1, 2, 2)
    plt.plot(time, width, color='b', marker= 'o')
    plt.title("Diffusion layer", fontsize=17)
    plt.grid()
    # plt.xlim(0, 0.4)
    plt.yscale('log')
    plt.ylabel('Width', fontsize=15, labelpad=5)
    plt.xlabel('time [s]', fontsize=15, labelpad=5)
    plt.tick_params(axis='both', which='major', labelsize=15)
    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)

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
d['mesh_1'+'_dirPath'] = '/home/aitor/Escritorio/Norouzi_mesh_1.gid/Norouzi_mesh_1_StationarityL2ErrorNorm.hdf5'
velocity_error_v, concentration_error_v, times_v = GetFieldDataFile(d['mesh_1'+'_dirPath'])
width_v = GetWidth(d['mesh_1'+'_dirPath'])
PlotMeshConvergence(velocity_error_v, concentration_error_v, times_v, width_v)

            




        