from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
import co_simulation_tools as cs_tools
from CoSimulationApplication import *
import sys
import from co_simulation_analysis import CoSimulationAnalysis

if __name__ == '__main__':
    from sys import argv

    if len(argv) != 2:
        err_msg =  'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]
    global cs_data_structure
    cs_data_structure = cs_tools.ImportDataStructure(parameter_file_name)

    # Now we import actual parameters from the cs_data_structure
    with open(parameter_file_name,'r') as parameter_file:
        parameters = cs_data_structure.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()

