from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.IgaApplication.iga_analysis import IgaAnalysis

"""
For user-scripting it is intended that a new class is derived
from IgaAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = IgaAnalysis(model,parameters)
    simulation.Run()
    
    # simulation.Initialize()

    # structural_analysis_model_part = model.GetModelPart("IgaModelPart.StructuralAnalysis")

    # for node in structural_analysis_model_part.Nodes:
    #     node.Fix(KratosMultiphysics.ROTATION_X)
    #     node.Fix(KratosMultiphysics.ROTATION_Y)
    #     node.Fix(KratosMultiphysics.ROTATION_Z)

    # simulation.RunSolutionLoop()
    # simulation.Finalize()

