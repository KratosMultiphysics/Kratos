from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
from KratosMultiphysics.mpi import *
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
KratosMultiphysics.CheckForPreviousImport()
import trilinos_navier_stokes_solver_vmsmonolithic
import mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return ALETrilinosNavierStokesSolverVMSMonolithic(model_part, custom_settings)


class ALETrilinosNavierStokesSolverVMSMonolithic(trilinos_navier_stokes_solver_vmsmonolithic.NavierStokesMPISolver_VMSMonolithic):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("Monolithic")
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).__init__(model_part, navier_stokes_settings)
        # create ale solver
        ale_solver_type = custom_settings["ale_settings"]["solver_type"].GetString()
        valid_ale_solver_types = ["trilinos_mesh_solver_structural_similarity"]
        if ale_solver_type not in valid_ale_solver_types:
            raise Exception("Invalid ALE solver_type: " + ale_solver_type)
        ale_solver_module = __import__(ale_solver_type)
        self.ale_solver = ale_solver_module.CreateSolver(self.main_model_part, custom_settings["ale_settings"])
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: Construction finished")

    def AddVariables(self):
        # add base class variables
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).AddVariables()
        self.ale_solver.AddVariables()
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: Variables ADDED.")

    def AddDofs(self):
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).AddDofs()
        self.ale_solver.AddDofs()
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: DOFs ADDED.")

    def Initialize(self):
        super(ALETrilinosNavierStokesSolverVMSMonolithic, self).Initialize()
        self.ale_solver.Initialize()
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverVMSMonolithic]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALETrilinosNavierStokesSolverVMSMonolithic, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        self.GetMeshMotionSolver().Solve()
        self.GetFluidSolver.Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
