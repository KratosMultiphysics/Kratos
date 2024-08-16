import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.inverseforming_analysis import InverseFormingAnalysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = InverseFormingAnalysis(model,parameters)
    simulation.Run()
