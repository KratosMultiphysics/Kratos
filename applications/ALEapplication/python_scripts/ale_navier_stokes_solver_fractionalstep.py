from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
KratosMultiphysics.CheckForPreviousImport()
import navier_stokes_solver_fractionalstep


def CreateSolver(model_part, custom_settings):
    return ALENavierStokesSolverFractionalStep(model_part, custom_settings)


class ALENavierStokesSolverFractionalStep(navier_stokes_solver_fractionalstep.NavierStokesSolver_FractionalStep):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        navier_stokes_settings = custom_settings.Clone()
        navier_stokes_settings.RemoveValue("ale_settings")
        navier_stokes_settings["solver_type"].SetString("FractionalStep")
        super(ALENavierStokesSolverFractionalStep, self).__init__(model_part, navier_stokes_settings)
        # create ale solver
        ale_solver_type = custom_settings["ale_settings"]["solver_type"].GetString()
        valid_ale_solver_types = ["mesh_solver_structural_similarity"]
        if ale_solver_type not in valid_ale_solver_types:
            raise Exception("Invalid ALE solver_type: " + ale_solver_type)
        ale_solver_module = __import__(ale_solver_type)
        self.ale_solver = ale_solver_module.CreateSolver(self.main_model_part, custom_settings["ale_settings"])
        print("::[ALENavierStokesSolverFractionalStep]:: Construction finished")

    def AddVariables(self):
        # add base class variables
        super(ALENavierStokesSolverFractionalStep, self).AddVariables()
        # add ale variables
        self.ale_solver.AddVariables()
        print("::[ALENavierStokesSolverFractionalStep]:: Variables ADDED")

    def AddDofs(self):
        # add base class dofs
        super(ALENavierStokesSolverFractionalStep, self).AddDofs()
        # add ale dofs
        self.ale_solver.AddDofs()
        print("::[ALENavierStokesSolverFractionalStep]:: DOFs ADDED")

    def Initialize(self):
        super(ALENavierStokesSolverFractionalStep, self).Initialize()
        self.ale_solver.Initialize()
        print("::[ALENavierStokesSolverFractionalStep]:: Finished initialization.")

    def GetFluidSolver(self):
        return super(ALENavierStokesSolverFractionalStep, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver

    def Solve(self):
        self.GetMeshMotionSolver().Solve()
        self.GetFluidSolver().Solve()

    def MoveMesh(self):
        self.GetMeshMotionSolver().MoveMesh()
