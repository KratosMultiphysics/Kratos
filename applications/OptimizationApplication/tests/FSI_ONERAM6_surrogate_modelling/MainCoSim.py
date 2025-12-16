
import KratosMultiphysics as KM
import time
#import KratosMultiphysics.EmpireApplication

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
import json

# parameter_file_name = "project_parameters_cosim_pure_SDoF.json"
# parameter_file_name = "project_parameters_cosim_pure_fluid.json"
parameter_file_name = "ProjectParametersFSI.json"

with open(parameter_file_name, 'r') as parameter_file:
    cosim_parameters = KM.Parameters(parameter_file.read())
start = time.time()
simulation = CoSimulationAnalysis(cosim_parameters)
simulation.Run()
end = time.time() - start
print("time:  ", end)

