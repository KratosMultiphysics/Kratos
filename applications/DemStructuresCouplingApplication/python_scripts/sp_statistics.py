import math
import matplotlib.pyplot as plt
import numpy as np
import h5py

# Read HDF5 file
file_name = '/home/ipouplana/Ara/DEM_FEM/ctw16_pdot5e10_friction65_cohesion0.0_1.5R/sp_data.hdf5'
f = h5py.File(file_name, 'r')

# File Attributes and general factors
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

height_factor = real_probe_height / thickness
porosity_factor = (1.0-target_porosity)/(1.0-porosity)

# Total initial mass
initial_mass = 0.0
gram_factor = 1000.0
initial_num_spheres = len(f['0/radius'])
# Spheres loop
for i in range(initial_num_spheres):
    initial_mass +=4.0/3.0 * np.pi * f['0/radius'][i] ** 3 * density * gram_factor * height_factor * porosity_factor

# SP
total_num_steps = len(f.keys())
failure_step = 175 # TODO
step = 5
num_steps = int(failure_step/step)

all_times=[]
all_pressures=[]
all_sps=[]

times = np.array([])
for i in range(0,failure_step,step):
    times = np.append(times,f[str(i)].attrs['time'])
pressures = [5e10*t for t in times] # Pressure in Pa
pressures = [0.000145038*p for p in pressures] # Pressure in psi

# times = np.zeros(num_steps)
# sps = np.zeros(num_steps)

# sps_0 = np.array([])
# sps_1 = np.array([])
# sps_2 = np.array([])
# sps_3 = np.array([])
# sps_4 = np.array([])
# sps_5 = np.array([])

for nb in np.arange(0.0,6.0,1.0):
    sps = np.array([])
    for i in range(0,failure_step,step):
        # times[int(i/num_steps)] = f[str(i)].attrs['time']
        times = np.append(times,f[str(i)].attrs['time'])
        current_mass = 0.0
        # TODO: use numpy arrays operations instead of standard loop (this is very slow!)
        # all_radius_h5 = f[str(i)].get('radius')
        # all_radius = np.array(all_radius_h5)
        # Spheres loop
        for j in range(len(f[str(i)+'/radius'])):
            if f[str(i)+'/current_continuum_bonds'][j] > nb:
                current_mass += 4.0/3.0 * np.pi * f[str(i)+'/radius'][j] ** 3 * density * gram_factor * height_factor * porosity_factor
        current_sp = initial_mass - current_mass
        # sps[int(i/num_steps)] = current_sp
        sps = np.append(sps,current_sp)
        print(i)
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
graph_name='sp_bonds.pdf'
graph_labels=['SP up to 0 intact bonds',
            'SP up to 1 intact bonds',
            'SP up to 2 intact bonds',
            'SP up to 3 intact bonds',
            'SP up to 4 intact bonds',
            'SP up to 5 intact bonds',
            'ctw16 experiment',
                ]

f = plt.figure()

# plt.plot(pressures, sps)
# plt.xlabel('Applied pressure (psi)')
# plt.ylabel('Sand Production (g)')
# plt.title('SP with 0 or 1 bond left')
# f.savefig("sp_0_1.pdf", bbox_inches='tight')

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