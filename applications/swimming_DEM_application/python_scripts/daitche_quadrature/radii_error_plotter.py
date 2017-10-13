import os
import matplotlib.pyplot as plt
import itertools
import math


def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

dir = os.path.dirname(__file__)
my_paths = []
my_labels = []
my_Nqs = []
my_file_path    = os.path.join(dir, 'Daitche_1_order_2Post_Files/radii1.txt')
my_paths += [my_file_path]
my_Nqs += [1]
my_labels += [r'$N_q = 1$']
my_file_path    = os.path.join(dir, 'Daitche_2_order_2Post_Files/radii2.txt')
#my_paths += [my_file_path]
#my_labels += [r'$N_q = 2$']
#my_Nqs += [2]
my_file_path    = os.path.join(dir, 'Daitche_4_order_2Post_Files/radii4.txt')
my_paths += [my_file_path]
my_labels += [r'$N_q = 4$']
my_Nqs += [4]
my_file_path    = os.path.join(dir, 'Daitche_8_order_2Post_Files/radii8.txt')
#my_paths += [my_file_path]
#my_labels += [r'$N_q = 8$']
#my_Nqs += [8]
my_file_path    = os.path.join(dir, 'Daitche_16_order_2Post_Files/radii16.txt')
my_paths += [my_file_path]
my_labels += [r'$N_q = 16$']
my_Nqs += [16]
my_file_path    = os.path.join(dir, 'Daitche_32_order_2Post_Files/radii32.txt')
#my_paths += [my_file_path]
#my_labels += [r'$N_q = 32$']
#my_Nqs += [32]
my_file_path    = os.path.join(dir, 'Daitche_64_order_2Post_Files/radii64.txt')
my_paths += [my_file_path]
my_labels += [r'$N_q = 64$']
my_Nqs += [64]
exact_file_path = os.path.join(dir, 'Reference_Post_Files/exact_radii.txt')
my_radii_arrays = []
my_labels 

for my_file_path in my_paths:
    with open(my_file_path, mode = 'r') as f:
        ts = []
        my_radii = []
        for line in f:
            data = str.split(line)
            ts.append(float(data[0]))
            my_radii.append(float(data[1]))
        my_radii_arrays += [my_radii]
        
with open(exact_file_path, mode = 'r') as f:
    exact_radii = []
    for line in f:
        data = str.split(line)
        exact_radii.append(float(data[1]))

exact_file_path_2 = os.path.join(dir, 'Reference_Post_Files_2/exact_radii.txt')
exact_file_path_4 = os.path.join(dir, 'Reference_Post_Files_4/exact_radii.txt')
exact_file_path_8 = os.path.join(dir, 'Reference_Post_Files_8/exact_radii.txt')
exact_file_path_16 = os.path.join(dir, 'Reference_Post_Files_16/exact_radii.txt')

with open(exact_file_path_2, mode = 'r') as f:
    exact_radii_2 = []    
    for line in f:
        data = str.split(line)
        exact_radii_2.append(float(data[1]))
        
with open(exact_file_path_4, mode = 'r') as f:
    exact_radii_4 = []    
    for line in f:
        data = str.split(line)
        exact_radii_4.append(float(data[1]))
        
with open(exact_file_path_8, mode = 'r') as f:
    exact_radii_8 = []    
    for line in f:
        data = str.split(line)
        exact_radii_8.append(float(data[1]))
        
with open(exact_file_path_16, mode = 'r') as f:
    exact_radii_16 = []    
    for line in f:
        data = str.split(line)
        exact_radii_16.append(float(data[1]))        

j = 0
for my_radii in my_radii_arrays:
    rel_error = [abs(exact_radii[i] - my_radii[i]) / exact_radii[i] for i in range(len(ts))]
    print(len(rel_error))
    plt.plot(ts, rel_error, label=my_labels[j], marker = 4 + j, markevery = int(1000), color='k')
    j += 1


my_file_path_2 = os.path.join(dir, 'Daitche_t_step_2Post_Files/radii1.txt')
with open(my_file_path_2, mode = 'r') as f:
    ts = []
    my_radii = []
    for line in f:
        data = str.split(line)
        ts.append(float(data[0]))
        my_radii.append(float(data[1]))
    my_radii_arrays += [my_radii]
rel_error = [abs(exact_radii_2[i] - my_radii[i]) / exact_radii_2[i] for i in range(len(ts))]
plt.plot(ts, rel_error, label=r'$\Delta t \times 2$', linestyle='--', marker = '*', markevery = int(1000 / 2), color='k')

my_file_path_4 = os.path.join(dir, 'Daitche_t_step_4Post_Files/radii1.txt')
with open(my_file_path_4, mode = 'r') as f:
    ts = []
    my_radii = []
    for line in f:
        data = str.split(line)
        ts.append(float(data[0]))
        my_radii.append(float(data[1]))
    my_radii_arrays += [my_radii]
rel_error = [abs(exact_radii_4[i] - my_radii[i]) / exact_radii_4[i] for i in range(len(ts))]
plt.plot(ts, rel_error, label=r'$\Delta t \times 4$', linestyle='--', marker = 5, markevery = int(1000 / 4), color='k')

my_file_path_8 = os.path.join(dir, 'Daitche_t_step_8Post_Files/radii1.txt')
with open(my_file_path_8, mode = 'r') as f:
    ts = []
    my_radii = []
    for line in f:
        data = str.split(line)
        ts.append(float(data[0]))
        my_radii.append(float(data[1]))
    my_radii_arrays += [my_radii]
rel_error = [abs(exact_radii_8[i] - my_radii[i]) / exact_radii_8[i] for i in range(len(ts))]
plt.plot(ts, rel_error, label=r'$\Delta t \times 8$', linestyle='--', marker = 3, markevery = int(1000 / 8), color='k')

my_file_path_16 = os.path.join(dir, 'Daitche_t_step_16Post_Files/radii1.txt')
with open(my_file_path_16, mode = 'r') as f:
    ts = []
    my_radii = []
    for line in f:
        data = str.split(line)
        ts.append(float(data[0]))
        my_radii.append(float(data[1]))
    my_radii_arrays += [my_radii]
rel_error = [abs(exact_radii_16[i] - my_radii[i]) / exact_radii_16[i] for i in range(len(ts))]
plt.plot(ts, rel_error, label=r'$\Delta t \times 16$', linestyle='--', marker = 6, markevery = int(1000 / 16), color='k')


plt.semilogy()
plt.legend(loc='bottom right')
plt.ylabel(r'$E_{\mathrm{rel}}$')
plt.xlabel(r'$t$')
ax = plt.subplot(111)
handles, labels = ax.get_legend_handles_labels()
plt.legend(flip(handles, 2), flip(labels, 2), loc=4, ncol=2)
plt.savefig('ErrorInTime.eps', format='eps', dpi=1200)
plt.savefig('ErrorInTime.pdf', format='pdf', dpi=1200)

plt.show()