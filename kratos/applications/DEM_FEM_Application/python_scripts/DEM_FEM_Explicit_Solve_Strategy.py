from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library


from KratosMultiphysics import *
from KratosMultiphysics.DEM_FEM_Application import *


def AddVariables(model_part):

    # Just for meshing, without any meanings
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(NODAL_H)

    print("Variables for the DEM&FEM Coupled solution added correctly")


def AddDofs(model_part):
    # for node in model_part.Nodes:

    print("Dofs for the DEM&FEM Coupled solution added correctly")


class ExplicitStrategy:
    #

    def __init__(self, fem_model_part, dem_model_part, domain_size, dem_strategy, fem_strategy):
        self.fem_model_part = fem_model_part
        self.dem_model_part = dem_model_part
        self.domain_size = domain_size
        self.damp_option = 1
        self.damp_ratio = 0.2
        self.MoveMeshFlag = True
        self.dem_strategy = dem_strategy
        self.fem_strategy = fem_strategy

    #
    def InitializeStrategy(self):

        # creating the solution strategy
        self.solver = DemFemExplicitSolverStrategy(self.fem_model_part, self.dem_model_part, self.domain_size, self.damp_option, self.damp_ratio, self.MoveMeshFlag, self.dem_strategy, self.fem_strategy)

    #
        def Solve(self):
            (self.solver).Solve()

        def Initialize(self):
            (self.solver).Initialize()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
