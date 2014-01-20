from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.DEM_FEM_Application import *
from KratosMultiphysics.StructuralApplication import *

CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(GRAVITY)
    model_part.AddNodalSolutionStepVariable(NODAL_STRESS)
    model_part.AddNodalSolutionStepVariable(VIRTUAL_NODAL_MASS);

    print("Variables for the dynamic structural solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(VELOCITY_X, REACTION_X);
        node.AddDof(VELOCITY_Y, REACTION_Y);
        node.AddDof(VELOCITY_Z, REACTION_Z);
    print("Dofs for the dynamic structural solution added correctly")


class DynamicFEMSolver:
    #

    def __init__(self, model_part, Param):

        self.model_part = model_part
        self.domain_size = int(Param.Dimension)
        self.damping_ratio = Param.FemDampRatio
        self.max_delta_time = Param.MaxTimeStep;
        self.CalculateReactionFlag = True;
        self.MoveMeshFlag = True;

        if(Param.FemDampType == "LocalDamp"):
            self.damp_type = 1
        else:
            self.damp_type = 2

        if(Param.VirtualMassOption == "ON"):
            self.If_Virtual_Mass = True
        else:
            self.If_Virtual_Mass = False

    #
    def Initialize(self):

        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.max_delta_time)

        # creating the solution strategy
        self.solver = FemExplicitForwardStrategy(self.model_part, self.domain_size, self.damp_type, self.damping_ratio, self.If_Virtual_Mass, self.max_delta_time, self.CalculateReactionFlag, self.MoveMeshFlag)

    #
    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
