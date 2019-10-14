from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

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

solver = cs_tools.CreateInstance(parameters["solver_wrappers"][0])

interface_input = solver.GetInterfaceInput()

solver.SolveSolutionStep(interface_input)
