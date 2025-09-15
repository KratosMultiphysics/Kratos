
import KratosMultiphysics

from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

"""
For user-scripting it is intended that a new class is derived
from PotentialFlowAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PotentialFlowAnalysis(model,parameters)
    simulation.Run()
