from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure

import numpy as np
from sys import argv
import os

def print_colored(string, color):
    if color=='green':
        print('\x1b[0;30;42m'+string+'\x1b[0m')
    elif color=='orange':
        print('\x1b[0;30;43m' + string + '\x1b[0m')
    elif color=='red':
        print('\x1b[0;30;41m' + string + '\x1b[0m')
    else:
        print(string+f'(color {color} not implemented)')


# Check number of command line arguments
if len(argv) != 2:
    err_msg = 'Wrong number of input arguments!\n'
    err_msg += 'Use this script in the following way:\n'
    err_msg += '    "python Test.py <parameter-file>.json"\n'
    raise Exception(err_msg)

# Import data structure
parameter_file_name = argv[1]
cs_data_structure = ImportDataStructure(parameter_file_name)

# Import parameters using the data structure
with open(parameter_file_name, 'r') as parameter_file:
    parameters = cs_data_structure.Parameters(parameter_file.read())

# Create the solver (__init__)
print("Creating an AbaqusSolver")
AbaqusSolver0  = cs_tools.CreateInstance(parameters["solver_wrappers"][0])
print_colored("AbaqusSolver0 created","green")

# Assign loads to the Input-Nodes
# give value to DISPLACEMENT variable
mp = AbaqusSolver0.model['BEAMINSIDEMOVING_load_points'] #interface input modelpart
pressure = vars(KM)['PRESSURE']
traction = vars(KM)['TRACTION']
p = 10000
for node in mp.Nodes:
    # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
    # print(node.Y)
    node.SetSolutionStepValue(pressure, 0, p)
    node.SetSolutionStepValue(traction, 0, [0,0,0])
print(f"Assigned uniform pressure ({p} Pa) and 0 traction at the interface ")

AbaqusSolver0.Initialize()

#Step 1, Coupling 1
AbaqusSolver0.InitializeSolutionStep()
AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())

os.system("cp -r CSM/CSM_Time1.odb CSM/CSM_Time1_Iter1.odb")

#Step 1, Coupling 2
AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
AbaqusSolver0.FinalizeSolutionStep()

#Step 2, Coupling 1
AbaqusSolver0.InitializeSolutionStep()
AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
AbaqusSolver0.FinalizeSolutionStep()

#Iterate until deformation is approximately steady
mp_out = AbaqusSolver0.model['BEAMINSIDEMOVING_nodes'] #interface input modelpart
displacement = vars(KM)['DISPLACEMENT']

n_out = mp_out.NumberOfNodes()
prev_displacement = np.zeros((n_out,3))*0.
diff = np.zeros((n_out,3))*0.
tol = 1e-06
diffMax = 1000
while diffMax > tol:
    AbaqusSolver0.InitializeSolutionStep()
    AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
    AbaqusSolver0.FinalizeSolutionStep()
    diffMax = 0
    for node in mp_out.Nodes:
        diff = np.linalg.norm(np.array(node.GetSolutionStepValue(displacement))-prev_displacement[int(node.Id),:])
        prev_displacement[int(node.Id),:] = np.array(node.GetSolutionStepValue(displacement))
        if diff > diffMax:
            diffMax = diff
    print(diffMax)


AbaqusSolver0.Finalize()

print_colored("Finished",'green')