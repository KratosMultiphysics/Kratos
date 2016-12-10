from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY_LAPLACIAN)

def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_LAPLACIAN_X)
        node.AddDof(VELOCITY_LAPLACIAN_Y)
        node.AddDof(VELOCITY_LAPLACIAN_Z)

    print("dofs for the derivative recovery solver added correctly")
