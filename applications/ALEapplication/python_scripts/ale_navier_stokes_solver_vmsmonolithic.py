from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ALEApplication", "FluidDynamicsApplication","ExternalSolversApplication")

# Import applications
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication

# Other imports
import os
import navier_stokes_solver_vmsmonolithic


def CreateSolver(model_part, custom_settings):
    return ALENavierStokesSolverVMSMonolithic(model_part, custom_settings)


class ALENavierStokesSolverVMSMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolver_VMSMonolithic):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("Monolithic")
        super(ALENavierStokesSolverVMSMonolithic, self).__init__(model_part, navier_stokes_settings)
        # create ale solver
        ale_solver_type = custom_settings["ale_settings"]["solver_type"].GetString()
        valid_ale_solver_types = ["mesh_solver_structural_similarity", "mesh_solver_laplacian"]
        if ale_solver_type not in valid_ale_solver_types:
            raise Exception("Invalid ALE solver_type: " + ale_solver_type)
        ale_solver_module = __import__(ale_solver_type)
        self.ale_solver = ale_solver_module.CreateSolver(self.main_model_part, custom_settings["ale_settings"])
        print("::[ALENavierStokesSolverVMSMonolithic]:: Construction finished")

    def AddVariables(self):
        super(ALENavierStokesSolverVMSMonolithic, self).AddVariables()
        self.ale_solver.AddVariables()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Variables ADDED.")

    def AddDofs(self):
        super(ALENavierStokesSolverVMSMonolithic, self).AddDofs()
        self.ale_solver.AddDofs()
        print("::[ALENavierStokesSolverVMSMonolithic]:: DOFs ADDED.")

    def Initialize(self):
        super(ALENavierStokesSolverVMSMonolithic, self).Initialize()
        self.ale_solver.Initialize()
        print("::[ALENavierStokesSolverVMSMonolithic]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALENavierStokesSolverVMSMonolithic, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        self.GetMeshMotionSolver().Solve()
        self.GetFluidSolver().Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
