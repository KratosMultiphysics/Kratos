from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.AdjointFluidApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(ADJOINT_VELOCITY)
    model_part.AddNodalSolutionStepVariable(ADJOINT_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(ADJOINT_PRESSURE)
    #model_part.AddNodalSolutionStepVariable(PRIMAL_VELOCITY)
    #model_part.AddNodalSolutionStepVariable(PRIMAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)
    model_part.AddNodalSolutionStepVariable(NORMAL_SENSITIVITY)
    print("variables for the adjoint fluid solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(ADJOINT_VELOCITY_X)
        node.AddDof(ADJOINT_VELOCITY_Y)
        node.AddDof(ADJOINT_VELOCITY_Z)
        node.AddDof(ADJOINT_PRESSURE)

    print("dofs for the adjoint fluid solver added correctly")


class AdjointFluidSolver:
    
    def __init__(self, model_part, dimension):

        self.model_part = model_part
        self.dimension = dimension

        import KratosMultiphysics.ExternalSolversApplication
        self.linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUSolver()
        #self.linear_solver = BICGSTABSolver(1e-9, 300)

        print("Construction adjoint fluid solver finished")

    def Initialize(self):
        self.solver = AdjointFluidStrategy(self.model_part, self.linear_solver, self.dimension)

        print ("Initialization adjoint fluid solver finished")
    
    def Solve(self):
        (self.solver).Solve()

    def SetDragForceDirection(self, direction):
        (self.solver).SetDragForceDirection(direction)

    def ComputeSensitivity(self):
        (self.solver).ComputeSensitivity()
