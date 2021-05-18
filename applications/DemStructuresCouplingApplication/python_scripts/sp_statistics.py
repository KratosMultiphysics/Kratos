import matplotlib.pyplot as plt
import numpy as np
import h5py

# Read HDF5 file
file_name = 'sp_data.hdf5'
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
density = f.attrs['density']

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
all_times=[]
all_pressures=[]
all_sps=[]

# Time and external applied pressure
failure_step = 175 # TODO
max_step = 2*failure_step
max_step = failure_step
times = np.zeros(max_step+1)
for i in range(max_step):
    times[i+1] = f[str(i)].attrs['time']
p_rate = 5e10 # TODO
psi_factor = 0.000145038
pressures = [p_rate*t*psi_factor for t in times] # Pressure in psi
t_f=times[failure_step]
# print(t_f)
# Estimated time of initial sanding. TODO
# t_is=times[90]
# t_is=times[130]
# print(t_is)

# Maximum number of intact bonds to consider spheres as SP
max_num_bonds = 6.0 # TODO: this should be an int

for numbonds in np.arange(0.0,max_num_bonds,1.0):
    sps = np.zeros(max_step+1)
    for i in range(max_step):

        # Read datasets
        all_radii = np.array(f[str(i)].get('radius'))
        continuum_bonds = np.array(f[str(i)].get('current_continuum_bonds'))
        xs = np.array(f[str(i)].get('x'))
        ys = np.array(f[str(i)].get('y'))

        # Separate internal and external spheres (weak region and strong region)
        xs_2 = np.power(xs,2)
        ys_2 = np.power(ys,2)
        distance_2 = xs_2 + ys_2
        # weak_radius = 0.5*(internal_radius+interface_radius)
        # if times[i+1] < t_is:
        #     weak_radius = internal_radius
        # else:
        #     weak_radius = internal_radius + (interface_radius-internal_radius)/(t_f-t_is)*(times[i+1]-t_is)
        # TODO: count all the spheres in the domain
        weak_radius = interface_radius
        internal_radii = np.where(distance_2<weak_radius**2,all_radii,0.0) # spheres between inner wall (hole) and a certain radius (weak region)
        external_radii = np.where(distance_2>=weak_radius**2,all_radii,0.0) # spheres between a certain radius and the interface wall (dem-fem wall) (strong region)

        # Eliminate spheres that are free (SP), taking special care of those falling near the hole (internal spheres)
        cont_internal_radii = np.where(continuum_bonds>numbonds,internal_radii,0.0) # We eliminate spheres with small number of bonds in the weak part of the probe
        cont_external_radii = np.where(continuum_bonds>0.0,external_radii,0.0)

        total_radii_3 = np.power(cont_internal_radii+cont_external_radii,3)

        # Compute SP
        masses = 4.0/3.0 * np.pi * density * gram_factor * height_factor * porosity_factor * total_radii_3
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
graph_name='sp_bonds_t.pdf'
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

for name, pressures, productions in zip(graph_labels, all_pressures, all_sps):
    plt.plot(pressures, productions,label=name)
#p_is = p_rate*t_is*psi_factor
#plt.axvline(x=p_is, c='k', ls='--', label='estimated numerical initial sanding')

# for name, times, productions in zip(graph_labels, all_times, all_sps):
#     plt.plot(times, productions,label=name)
# plt.axvline(x=t_is, c='k', ls='--', label='estimated numerical initial sanding')
# plt.axvline(x=t_f, c='k', ls=':', label='numerical collapse')

plt.legend(loc=2, prop={'size': 6})
plt.xlabel('p (psi)')
# plt.xlabel('Time (s)')
plt.ylabel('Sand Production (g)')
plt.title('SP depending on number of intact bonds')
f.savefig(graph_name, bbox_inches='tight')