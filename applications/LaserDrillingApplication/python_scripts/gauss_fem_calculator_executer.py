import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import random
from fem_tools import SurfaceFEMProjector
import time as timer
from os import environ
environ['OMP_NUM_THREADS'] = '4'

starting_time = timer.time()
n_elements = 10 # number of FEM elements
R_far = 1.0
debug_mode = True
H = 20
sparse_option = True
number_of_triangles = 110 # number of material elements
n_evap_elements = number_of_triangles
random.seed(42)

#delta_coefficients = np.array([0.5 * R_far * ((i+1)%2) * random.random() * i / n_evap_elements for i in range(2 * n_evap_elements)])
#evap_element_centers = np.array([0.5 * R_far * ((i+1)%2) * random.random() * i / n_evap_elements for i in range(2 * n_evap_elements)])
#evap_enthalpies = np.array([r**2 for r in evap_element_centers])
#evap_element_centers = np.array([0.25 * R_far * (i+1) / n_evap_elements for i in range(4 * n_evap_elements)])

evap_element_centers = np.array([R_far * random.random() for i in range(n_evap_elements)])
evap_element_centers.sort()

elems_sides = R_far / len(evap_element_centers)

evap_enthalpies = np.array([H * elems_sides * r for r in evap_element_centers])
support_elements = [[] for i in range(n_elements+1)]

projector = SurfaceFEMProjector(n_elements, R_far, sparse_option) #, delta_coefficients)

if not sparse_option:
    projector.FillUpMassMatrix()
else:
    projector.FillUpSparseMassMatrix()

projector.AssignDeltasToTestFunctionSupports(evap_element_centers, support_elements)
projector.FillUpDeltasRHS(evap_element_centers, support_elements, evap_enthalpies)
u = projector.Project()
q_cont = np.array([projector.q(r) for r in projector.X])
q_interp = projector.InterpolateFunctionAndNormalize(projector.q) #, 1.0)
remaining_energy = q_interp - u
q_eval = np.array([projector.EvaluateFEMFunction(remaining_energy, x) for x in evap_element_centers])

total_energy = projector.CalculateEnergyOfFEMFunction(u)
u_mean = np.mean(u)

if debug_mode:
    #print('A =', projector.A)
    #print('b =', projector.b)
    #print('r-coordinates =', evap_element_centers)
    print('total energy exp. =', sum(evap_enthalpies))
    print('total energy calc. =', total_energy)
    print('u_mean =', u_mean)
    #print('radii: ', evap_element_centers)    
    #print('volumes: ', evap_enthalpies)
    #print('u: ', u)
    #print('check evaluation:', projector.EvaluateFEMFunction([2 for i in range(n_elements+1)], 1*R_far))
    '''def N_test(x):
        return projector.N(10, x)
    u_interp = projector.InterpolateFunctionAndNormalize(N_test, 1.0)'''

plt.grid()
plt.plot(projector.X, u, color='blue', marker='o')
#plt.plot(projector.X, q_cont, color='black', marker='x')
plt.plot(projector.X, q_interp, color='red', marker='+')
plt.plot(projector.X, remaining_energy, color='green', marker='*')
#plt.plot(evap_element_centers, q_eval, color='black', marker='x')

plt.legend(["fluence (lost)", "fluence (interpolated)", "fluence (remaining)"], loc="upper right")

elapsed_time = timer.time() - starting_time
print('\nTime: ', elapsed_time, 's\n')

plt.show()
