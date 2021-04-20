from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import bisect as bi
import math
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import os
import random

def GetFieldDataFile(mdpa_name):
    model_part_name = mdpa_name
    hf = h5py.File(model_part_name, 'r')
    velocity_error_vector = []
    concentration_error_vector = []
    group_names = list(hf.keys())
    group_names = list(map(int,group_names))
    group_names.sort(reverse=False)
    for i in (group_names):
        time = float(hf['/' + str(i)].attrs['time'])
        n = 2
        n_round = (math.floor(time * 10 ** n) / 10 ** n)
        n_elements = int(hf['/' + str(i)].attrs['n_elements'])
        n_size = float(hf['/' + str(i)].attrs['element_size'])
        if n_round == 0.4:
            velocity_error_vector = np.append(velocity_error_vector, hf[str(group_names[i-1])+'/VELOCITY_ERROR'])
            concentration_error_vector = np.append(concentration_error_vector, hf[str(group_names[i-1])+'/CONCENTRATION_ERROR'])

    return velocity_error_vector, concentration_error_vector, n_elements, n_size

def PlotMeshConvergence(velocity, concentration, n_elem, n_size):
            
    # plt.subplot(1, 2, 1)
    plt.plot(n_size[::-1], velocity[::-1], color='b', marker = 'o')
    #convergence rate
    # q1 = (math.log(velocity[0][0]/velocity[1][0])/math.log(n_size[0]/n_size[1]))
    # q0 = (math.log(velocity[1][0]/velocity[2][0])/math.log(n_size[1]/n_size[2]))
    plt.title("Velocity", fontsize=17)
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

    # plt.subplot(1, 2, 2)
    # plt.plot(n_elem[0], concentration[0], color='b', marker='o', label = '$Concentration$ mesh $2$')
    # plt.plot(n_elem[1], concentration[1], color='k', marker='+', label = '$Concentration$ mesh $3$')
    # plt.plot(n_elem[2], concentration[2], color='g', marker='^', label = '$Concentration$ mesh $4$')
    # plt.title("Concentration", fontsize=17)
    # plt.grid()
    # # plt.xlim(0, 11)
    # plt.yscale('log')
    # # plt.ylabel('L2  error  norm', fontsize=20, labelpad=25)
    # plt.xlabel('steps', fontsize=15, labelpad=5)
    # lgnd = plt.legend(loc = 'lower left', prop={'size':15}, frameon=False)
    # plt.tick_params(axis='both', which='major', labelsize=15)
    # # plt.semilogy()
    # ax = plt.gca()
    # ax.tick_params(axis='x', pad=20)
    # ax.tick_params(axis='y', pad=10)
            
    plt.suptitle("Mesh Convergence at 4/7 $T_res$", fontsize=20)
    figure = plt.gcf() # get current figure
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
    figure.set_size_inches(14, 11)
    plt.savefig('L2 error norm' + '.pdf', format='pdf', dpi=1000)
    plt.show()

d = {}
n_element = []
n_size = []
velocity = []
concentration = []
for j in range(2,5):
    d['mesh_'+str(j)+'_dirPath'] = '/home/aitor/Escritorio/Norouzi_mesh_'+str(j)+'.gid/Norouzi_mesh_'+str(j)+'.hdf5'
    velocity_error_s, concentration_error_s, n_element_s, n_size_s = GetFieldDataFile(d['mesh_'+str(j)+'_dirPath'])
    velocity.append(velocity_error_s)
    concentration.append(concentration_error_s)
    n_size.append(n_size_s)
    n_element.append(n_element_s)
    # ad0 = (n_element[j-2]*n_size[j-2]**2)*0.5

PlotMeshConvergence(velocity, concentration, n_element, n_size)

        
