import math
import matplotlib.pyplot as plt
import numpy as np
import h5py

# Read HDF5 file
file_name = '/home/ipouplana/Ara/DEM_FEM/ctw16_pdot5e10_friction65_cohesion0.0_1.5R/sp_data.hdf5'
f = h5py.File(file_name, 'r')

# File Attributes
test_id = f.attrs['test_id']
internal_radius = f.attrs['internal_radius']
external_radius = f.attrs['external_radius']
interface_radius = f.attrs['interface_radius']
thickness = f.attrs['thickness']
volume = f.attrs['volume']
real_probe_height = f.attrs['real_probe_height']
target_porosity = f.attrs['target_porosity']
porosity = f.attrs['porosity']
# density = f.attrs['density']
density = 2575.0 # TODO: provisional

# General factors
height_factor = real_probe_height / thickness
porosity_factor = (1.0-target_porosity)/(1.0-porosity)
gram_factor = 1000.0

# Total initial mass
initial_radii = np.array(f['0'].get('radius'))
initial_radii_3 = np.power(initial_radii,3)
initial_masses = 4.0/3.0 * np.pi * density * gram_factor * height_factor * porosity_factor * initial_radii_3
initial_mass = np.sum(initial_masses)
# print(initial_mass)

# SP
failure_step = 175 # TODO

times = np.zeros(failure_step+1)
for i in range(failure_step):
    times[i+1] = f[str(i)].attrs['time']
pressures = [5e10*t*0.000145038 for t in times] # Pressure in psi
# print(times)
# print(pressures)

total_num_bonds = 6.0 # TODO: this should be an int

all_times=[]
all_pressures=[]
all_sps=[]

for numbonds in np.arange(0.0,total_num_bonds,1.0):
    sps = np.zeros(failure_step+1)
    for i in range(failure_step):
        radii = np.array(f[str(i)].get('radius'))
        continuum_bonds = np.array(f[str(i)].get('current_continuum_bonds'))
        r_cb = np.vstack((radii,continuum_bonds))
        radii_3 = np.where(r_cb[1,:]>numbonds,np.power(r_cb[0,:],3),0.0)
        masses = 4.0/3.0 * np.pi * density * gram_factor * height_factor * porosity_factor * radii_3
        current_mass = np.sum(masses)
        current_sp = initial_mass - current_mass
        sps[i+1] = current_sp
    all_times.append(times)
    all_pressures.append(pressures)
    all_sps.append(sps)

# Experiment data
with open('ctw16_experiment.txt') as exp_f:
    exp_times=[]
    exp_pressures=[]
    exp_sps=[]
    for line in exp_f:
        fields = line.strip().split()
        if fields:
            exp_times.append(float(fields[0]))
            exp_pressures.append(float(fields[1]))
            exp_sps.append(float(fields[2]))
all_times.append(exp_times)
all_pressures.append(exp_pressures)
all_sps.append(exp_sps)

# Graphs
graph_name='sp_bonds_b.pdf'
graph_labels=['SP up to 0 intact bonds',
            'SP up to 1 intact bonds',
            'SP up to 2 intact bonds',
            'SP up to 3 intact bonds',
            'SP up to 4 intact bonds',
            'SP up to 5 intact bonds',
            'ctw16 experiment',
                ]

f = plt.figure()

for name, pressures, productions in zip(graph_labels, all_pressures, all_sps):
    plt.plot(pressures, productions,label=name)
# for name, pressures, productions, i_collapse in zip(graph_labels, all_files_pressures, all_files_productions, all_files_i_collapse):
#     plt.plot(pressures[:i_collapse], productions[:i_collapse],label=name)
# for name, times, productions in zip(graph_labels, all_files_times, all_files_productions):
#     plt.plot(times, productions,label=name)

plt.legend(loc=2, prop={'size': 6})
plt.xlabel('p (psi)')
# plt.xlabel('Time (s)')
plt.ylabel('Sand Production (g)')
plt.title('SP depending on number of intact bonds')
f.savefig(graph_name, bbox_inches='tight')