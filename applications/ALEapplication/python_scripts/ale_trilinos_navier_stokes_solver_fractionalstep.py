from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
from KratosMultiphysics.mpi import *
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
KratosMultiphysics.CheckForPreviousImport()
import trilinos_navier_stokes_solver_fractionalstep
import mesh_solver_base

def CreateSolver(model_part, custom_settings):
    return ALETrilinosNavierStokesSolverFractionalStep(model_part, custom_settings)


class ALETrilinosNavierStokesSolverFractionalStep(trilinos_navier_stokes_solver_fractionalstep.Trilinos_NavierStokesSolver_FractionalStep):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("FractionalStep")
        super(ALETrilinosNavierStokesSolverFractionalStep, self).__init__(model_part, navier_stokes_settings)
        # create ale solver
        ale_solver_type = custom_settings["ale_solver_type"].GetString()
        valid_ale_solver_types = ["trilinos_mesh_solver_structural_similarity"]
        if self.ale_solver_type not in valid_ale_solver_types:
            raise Exception("Invalid ale_solver_type: " + self.ale_solver_type)

        ale_solver_module = __import__(ale_solver_type)
        self.ale_solver = ale_solver_module.CreateSolver(self.main_model_part, custom_settings["ale_settings"])
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverFractionalStep]:: Construction finished")

    def AddVariables(self):
        # add base class variables
        super(ALETrilinosNavierStokesSolverFractionalStep, self).AddVariables()
        # add ale variables
        self.ale_solver.AddVariables()
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverFractionalStep]:: Variables ADDED.")

    def AddDofs(self):
        # add base class dofs
        super(ALETrilinosNavierStokesSolverFractionalStep, self).AddDofs()
        # add ale dofs
        self.ale_solver.AddDofs()
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverFractionalStep]:: DOFs ADDED.")

    def Initialize(self):
        super(ALETrilinosNavierStokesSolverFractionalStep, self).Initialize()
        self.ale_solver.Initialize()
        if mpi.rank == 0:
            print("::[ALETrilinosNavierStokesSolverFractionalStep]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALETrilinosNavierStokesSolverFractionalStep, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        self.GetMeshMotionSolver().Solve()
        self.GetFluidSolver().Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
