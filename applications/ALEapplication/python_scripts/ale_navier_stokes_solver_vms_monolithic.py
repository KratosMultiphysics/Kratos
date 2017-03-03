from __future__ import print_function, absolute_import, division
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
CheckForPreviousImport()

import navier_stokes_solver_vmsmonolithic

def CreateSolver(model_part, custom_settings):
    return ALENavierStokesSolverVMSMonolithic(model_part, custom_settings)

class ALENavierStokesSolverVMSMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolver_VMSMonolithic):

    def __init__(self, model_part, custom_settings):
        # get ale solver type
        self.ale_solver_type = custom_settings["ale_solver_type"].GetString()
        # remove the ale solver type from settings so we can reuse the navier stokes constructor
        settings = custom_settings.PrettyPrintJsonString()
        i_start = settings[:settings.find('"ale_solver_type"')].rfind("\n")
        i_end = settings.find("\n",i_start+1)
        settings = settings[:i_start] + settings[i_end:]
        # call navier stokes constructor
        super().__init__(model_part, custom_settings)
        # create ale solver
        ale_solver_module = __import__(self.ale_solver_type)
        ale_settings = Parameters("") # for now use default settings
        self.ale_solver = ale_solver_module.CreateSolver(self.main_model_part, ale_settings)
        print("Construction of ALENavierStokesSolverVMSMonolithic finished")

    def AddVariables(self):
        # add base class variables
        super().AddVariables()
        # add ale variables
        self.ale_solver.AddVariables()
        print("ALE variables added correctly")

    def AddDofs(self):
        # add base class dofs
        super().AddDofs()
        # add ale dofs
        self.ale_solver.AddDofs()
        print("ALE dofs added correctly")

    def Initialize(self):
        super().Initialize()
        self.ale_solver.Initialize()

    def SolveFluid(self):
        super().Solve()

    def SolveMeshMotion(self):
        self.ale_solver.Solve()

    def Solve(self):
        self.SolveMeshMotion()
        self.SolveFluid()

    def MoveNodes(self):
        self.ale_solver.MoveNodes()
