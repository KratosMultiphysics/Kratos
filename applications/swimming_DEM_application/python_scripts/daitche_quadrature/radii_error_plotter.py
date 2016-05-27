import os
import matplotlib.pyplot as plt
import math
dir = os.path.dirname(__file__)
my_file_path_1 = os.path.join(dir, 'Daitche_1_Post_Files/radii1.txt')
my_file_path_1000 = os.path.join(dir, 'Daitche_1000_Post_Files/radii1000.txt')
exact_file_path = os.path.join(dir, 'Reference_Post_Files/exact_radii.txt')
ts = []
my_radii_1 = []
my_radii_1000 = []
exact_radii = []

with open(my_file_path_1, mode = 'r') as f:
    for line in f:
        data = str.split(line)
        ts.append(float(data[0]))
        my_radii_1.append(float(data[1]))
        
with open(my_file_path_1000, mode = 'r') as f:
    for line in f:
        data = str.split(line)
        my_radii_1000.append(float(data[1]))

with open(exact_file_path, mode = 'r') as f:
    for line in f:
        data = str.split(line)
        exact_radii.append(float(data[1]))

rel_error_1 = [abs(exact_radii[i] - my_radii_1[i]) / exact_radii[i] for i in range(len(ts))]
rel_error_1000 = [abs(exact_radii[i] - my_radii_1000[i]) / exact_radii[i] for i in range(len(ts))]
plt.plot(ts, rel_error_1)
plt.plot(ts, rel_error_1000)
plt.semilogy()
plt.show()