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

# Time and external applied pressure
failure_step = 300 # TODO
max_step = failure_step
times = np.zeros(max_step)
for i in range(max_step):
    times[i] = f[str(i)].attrs['time']

delta_probe_radius = (interface_radius-internal_radius)/5.0
# internal_radius_2 = internal_radius**2
probe_radius_1 = (internal_radius+delta_probe_radius*1)**2
probe_radius_2 = (internal_radius+delta_probe_radius*2)**2
probe_radius_3 = (internal_radius+delta_probe_radius*3)**2
probe_radius_4 = (internal_radius+delta_probe_radius*4)**2
interface_radius_2 = interface_radius**2
avg_num_intact_bonds_0 = np.zeros(max_step)
avg_num_intact_bonds_1 = np.zeros(max_step)
avg_num_intact_bonds_2 = np.zeros(max_step)
avg_num_intact_bonds_3 = np.zeros(max_step)
avg_num_intact_bonds_4 = np.zeros(max_step)

avg_num_broken_bonds_0 = np.zeros(max_step)
avg_num_broken_bonds_1 = np.zeros(max_step)
avg_num_broken_bonds_2 = np.zeros(max_step)
avg_num_broken_bonds_3 = np.zeros(max_step)
avg_num_broken_bonds_4 = np.zeros(max_step)

for i in range(max_step):

    # Read datasets
    continuum_bonds = np.array(f[str(i)].get('current_continuum_bonds'))
    initial_continuum_bonds = np.array(f[str(i)].get('initial_continuum_bonds'))
    xs = np.array(f[str(i)].get('x'))
    ys = np.array(f[str(i)].get('y'))

    xs_2 = np.power(xs,2)
    ys_2 = np.power(ys,2)
    distance_2 = xs_2 + ys_2

    num_intact_bonds_0 = np.where(distance_2<probe_radius_1,continuum_bonds,-100.0)
    index_to_delete = np.where(num_intact_bonds_0<0.0)[0]
    num_intact_bonds_0 = np.delete(num_intact_bonds_0,index_to_delete)
    num_init_bonds_0 = np.where(distance_2<probe_radius_1,initial_continuum_bonds,-100.0)
    index_to_delete = np.where(num_init_bonds_0<0.0)[0]
    num_init_bonds_0 = np.delete(num_init_bonds_0,index_to_delete)
    num_broken_bonds_0 = num_init_bonds_0 - num_intact_bonds_0
    avg_num_intact_bonds_0[i] = np.sum(num_intact_bonds_0)/len(num_intact_bonds_0)
    avg_num_broken_bonds_0[i] = np.sum(num_broken_bonds_0)/len(num_broken_bonds_0)

    num_intact_bonds_1 = np.where(distance_2>=probe_radius_1,continuum_bonds,-100)
    num_intact_bonds_1 = np.where(distance_2<probe_radius_2,num_intact_bonds_1,-100)
    index_to_delete = np.where(num_intact_bonds_1<0.0)[0]
    num_intact_bonds_1 = np.delete(num_intact_bonds_1,index_to_delete)
    num_init_bonds_1 = np.where(distance_2>=probe_radius_1,initial_continuum_bonds,-100)
    num_init_bonds_1 = np.where(distance_2<probe_radius_2,num_init_bonds_1,-100)
    index_to_delete = np.where(num_init_bonds_1<0.0)[0]
    num_init_bonds_1 = np.delete(num_init_bonds_1,index_to_delete)
    num_broken_bonds_1 = num_init_bonds_1 - num_intact_bonds_1
    avg_num_intact_bonds_1[i] = np.sum(num_intact_bonds_1)/len(num_intact_bonds_1)
    avg_num_broken_bonds_1[i] = np.sum(num_broken_bonds_1)/len(num_broken_bonds_1)

    num_intact_bonds_2 = np.where(distance_2>=probe_radius_2,continuum_bonds,-100)
    num_intact_bonds_2 = np.where(distance_2<probe_radius_3,num_intact_bonds_2,-100)
    index_to_delete = np.where(num_intact_bonds_2<0.0)[0]
    num_intact_bonds_2 = np.delete(num_intact_bonds_2,index_to_delete)
    num_init_bonds_2 = np.where(distance_2>=probe_radius_2,initial_continuum_bonds,-100)
    num_init_bonds_2 = np.where(distance_2<probe_radius_3,num_init_bonds_2,-100)
    index_to_delete = np.where(num_init_bonds_2<0.0)[0]
    num_init_bonds_2 = np.delete(num_init_bonds_2,index_to_delete)
    num_broken_bonds_2 = num_init_bonds_2 - num_intact_bonds_2
    avg_num_intact_bonds_2[i] = np.sum(num_intact_bonds_2)/len(num_intact_bonds_2)
    avg_num_broken_bonds_2[i] = np.sum(num_broken_bonds_2)/len(num_broken_bonds_2)

    num_intact_bonds_3 = np.where(distance_2>=probe_radius_3,continuum_bonds,-100)
    num_intact_bonds_3 = np.where(distance_2<probe_radius_4,num_intact_bonds_3,-100)
    index_to_delete = np.where(num_intact_bonds_3<0.0)[0]
    num_intact_bonds_3 = np.delete(num_intact_bonds_3,index_to_delete)
    num_init_bonds_3 = np.where(distance_2>=probe_radius_3,initial_continuum_bonds,-100)
    num_init_bonds_3 = np.where(distance_2<probe_radius_4,num_init_bonds_3,-100)
    index_to_delete = np.where(num_init_bonds_3<0.0)[0]
    num_init_bonds_3 = np.delete(num_init_bonds_3,index_to_delete)
    num_broken_bonds_3 = num_init_bonds_3 - num_intact_bonds_3
    avg_num_intact_bonds_3[i] = np.sum(num_intact_bonds_3)/len(num_intact_bonds_3)
    avg_num_broken_bonds_3[i] = np.sum(num_broken_bonds_3)/len(num_broken_bonds_3)

    num_intact_bonds_4 = np.where(distance_2>probe_radius_4,continuum_bonds,-100.0)
    index_to_delete = np.where(num_intact_bonds_4<0.0)[0]
    num_intact_bonds_4 = np.delete(num_intact_bonds_4,index_to_delete)
    num_init_bonds_4 = np.where(distance_2>probe_radius_4,initial_continuum_bonds,-100.0)
    index_to_delete = np.where(num_init_bonds_4<0.0)[0]
    num_init_bonds_4 = np.delete(num_init_bonds_4,index_to_delete)
    num_broken_bonds_4 = num_init_bonds_4 - num_intact_bonds_4
    avg_num_intact_bonds_4[i] = np.sum(num_intact_bonds_4)/len(num_intact_bonds_4)
    avg_num_broken_bonds_4[i] = np.sum(num_broken_bonds_4)/len(num_broken_bonds_4)


# Graphs
fig1, axs1 = plt.subplots(3, 2, figsize=(7, 7))

graph_name='average_number_intact_bonds_probe_radius.pdf'

axs1[0,0].plot(times,avg_num_intact_bonds_0)
axs1[1,0].plot(times,avg_num_intact_bonds_1)
axs1[2,0].plot(times,avg_num_intact_bonds_2)
axs1[0,1].plot(times,avg_num_intact_bonds_3)
axs1[1,1].plot(times,avg_num_intact_bonds_4)

fig1.savefig(graph_name)

fig2, axs2 = plt.subplots(3, 2, figsize=(7, 7))

graph_name='average_number_broken_bonds_probe_radius.pdf'

axs2[0,0].plot(times,avg_num_broken_bonds_0)
axs2[1,0].plot(times,avg_num_broken_bonds_1)
axs2[2,0].plot(times,avg_num_broken_bonds_2)
axs2[0,1].plot(times,avg_num_broken_bonds_3)
axs2[1,1].plot(times,avg_num_broken_bonds_4)

fig2.savefig(graph_name)