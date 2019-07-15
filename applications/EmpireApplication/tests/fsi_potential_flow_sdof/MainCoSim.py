from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.EmpireApplication.co_simulation_steady_analysis import CoSimulationSteadyAnalysis
import json

parameter_file_name = "fsi_potential_flow_sdof/project_cosim_naca0012_small_fsi_parameters.json"

with open(parameter_file_name, 'r') as parameter_file:
    cosim_parameters = json.load(parameter_file)

simulation = CoSimulationSteadyAnalysis(cosim_parameters)
simulation.Run()
