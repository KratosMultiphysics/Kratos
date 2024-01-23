
import KratosMultiphysics
from KratosMultiphysics.MPMApplication.mpm_analysis import MPMAnalysis

"""
For user-scripting it is intended that a new class is derived
from MPMAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = MPMAnalysis(model,parameters)
    simulation.Run()
