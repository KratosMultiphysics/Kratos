# # Import system python
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math


def GetWidth(mdpa_name):
    model_part_name = mdpa_name
    hf = h5py.File(model_part_name, 'r')
    target_coordinate_x = 0.024
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
        time = float(hf['/' + str(i)].attrs['time'])
        n = 2
        n_round = (math.floor(time * 10 ** n) / 10 ** n)
        if n_round == 0.001:
            n_size = float(hf['/' + str(i)].attrs['element_size'])
            coordinates = hf[str(group_names[i-1])+'/COORDINATES']
            nodes_id = hf[str(group_names[i-1])+'/ID']
            # velocity = hf[str(group_names[i-1])+'/VELOCITIES']
            concentration = hf[str(group_names[i-1])+'/TEMPERATURES']
            id_target = np.where(coordinates[:,0] == target_coordinate_x)
            target_tolerance = [coordinates[0][id_target]-n_size, coordinates[0][id_target]+n_size]
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
    return width, n_size

def PlotWidth(width, n_size):
            
    # plt.subplot(1, 2, 1)
    plt.plot(n_size[::-1], width[::-1], color='b', marker = 'o')
    plt.title("Width", fontsize=17)
    plt.grid()
    plt.xlim(15e-5, 4e-5)
    plt.ylabel('L2  error  norm', fontsize=15, labelpad=5)
    plt.xlabel('mesh size', fontsize=15, labelpad=5)
    lgnd = plt.legend(loc = 'lower left', prop={'size':15}, frameon=False)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.semilogy()
    ax = plt.gca()
    ax.tick_params(axis='x', pad=20)
    ax.tick_params(axis='y', pad=10)
    
    plt.suptitle("W(h)", fontsize=20)
    figure = plt.gcf() # get current figure
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
    figure.set_size_inches(14, 11)
    plt.savefig('Width_h' + '.pdf', format='pdf', dpi=1000)
    plt.show()

d = {}
n_size = []
width = []
for j in range(1,5):
    d['mesh_'+str(j)+'_dirPath'] = '/home/aitor/Escritorio/Norouzi_mesh_'+str(j)+'.gid/Norouzi_mesh_'+str(j)+'.hdf5'
    # ad0 = (n_element[j-2]*n_size[j-2]**2)*0.5
    width_s, n_size_s = GetWidth(d['mesh_'+str(j)+'_dirPath'])
    n_size.append(n_size_s)
    width.append(width_s)

PlotWidth(width, n_size)
