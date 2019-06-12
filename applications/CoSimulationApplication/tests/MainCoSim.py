from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

parameter_file_name = "fsi_mok/cosim_mok_fsi_parameters.json"
global cs_data_structure
cs_data_structure = cs_tools.ImportDataStructure(parameter_file_name)

# Now we import actual parameters from the cs_data_structure
with open(parameter_file_name,'r') as parameter_file:
    # import json
    # parameters = json.load(parameter_file)
    parameters = cs_data_structure.Parameters(parameter_file.read())

model = cs_data_structure.Model()

simulation = CoSimulationAnalysis(model, parameters)
simulation.Run()