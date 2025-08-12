import sys
import time

import KratosMultiphysics as KM
import KratosMultiphysics.DEMApplication 
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
"""
For user-scripting it is intended that a new class is derived
from ParticleMechanicsAnalysis to do modifications
"""

class my_analysis(CoSimulationAnalysis):
    def OutputSolutionStep(self):
        super().OutputSolutionStep()

parameter_file_name = "compression_dem_fem.json"
with open(parameter_file_name,'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

simulation = my_analysis(parameters)
simulation.Run()


