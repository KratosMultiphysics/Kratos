from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np

from sys import argv

# Check number of command line arguments
if len(argv) != 2:
    err_msg = 'Wrong number of input arguments!\n'
    err_msg += 'Use this script in the following way:\n'
    err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
    raise Exception(err_msg)


# Import data structure
parameter_file_name = argv[1]
cs_data_structure = ImportDataStructure(parameter_file_name)

# Import parameters using the data structure
with open(parameter_file_name, 'r') as parameter_file:
    parameters = cs_data_structure.Parameters(parameter_file.read())

solver = cs_tools.CreateInstance(parameters['solver_wrappers'][0])


settings = parameters['solver_wrappers'][0]['settings']

# steady test
if 0:
    solver.Initialize()
    solver.InitializeSolutionStep()

    interface_input = solver.GetInterfaceInput()
    for iteration in range(3):
        iteration += 1
        print(f'\niteration {iteration}')
        solver.SolveSolutionStep(interface_input)
        interface_input = solver.GetInterfaceInput()
        for key in settings['interface_input'].keys():
            for node in interface_input.model[key].Nodes:
                dy = (1 - np.cos(2 * np.pi * node.X)) * 0.5 * 0.01
                node.SetSolutionStepValue(vars(KM)['DISPLACEMENT'], 0, [0., dy, 0.])

    solver.FinalizeSolutionStep()
    solver.Finalize()

# unsteady test
else:
    solver.Initialize()

    interface_input = solver.GetInterfaceInput()
    for timestep in range(1, 5):
        f = 0.005 * (-1) ** (timestep + 1)
        f = 0.05
        solver.InitializeSolutionStep()
        for iteration in range(1, 5):
            solver.SolveSolutionStep(interface_input)
            interface_input = solver.GetInterfaceInput()
            for key in settings['interface_input'].keys():
                for node in interface_input.model[key].Nodes:
                    dy = (1 - np.cos(2 * np.pi * (node.X - timestep / 4 - iteration / 16))) * 0.5 * f
                    node.SetSolutionStepValue(vars(KM)['DISPLACEMENT'], 0, [0., dy, 0.])
        solver.FinalizeSolutionStep()

    solver.Finalize()