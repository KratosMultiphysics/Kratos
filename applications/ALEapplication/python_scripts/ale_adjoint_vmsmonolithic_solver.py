from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
KratosMultiphysics.CheckForPreviousImport()
import adjoint_vmsmonolithic_solver


def CreateSolver(model_part, custom_settings):
    return ALEAdjointVMSMonolithicSolver(model_part, custom_settings)


class ALEAdjointVMSMonolithicSolver(adjoint_vmsmonolithic_solver.AdjointVMSMonolithicSolver):
    def __init__(self, model_part, custom_settings):
        # remove the ale_settings so we can use the navier stokes constructor
        adjoint_settings = custom_settings.Clone()
        adjoint_settings.RemoveValue("ale_settings")
        super(ALEAdjointVMSMonolithicSolver, self).__init__(model_part, adjoint_settings)
        # create ale solver
        ale_solver_type = custom_settings["ale_settings"]["solver_type"].GetString()
        valid_ale_solver_types = ["mesh_solver_structural_similarity"]
        if ale_solver_type not in valid_ale_solver_types:
            raise Exception("Invalid ALE solver_type: " + ale_solver_type)
        ale_solver_module = __import__(ale_solver_type)
        self.ale_solver = ale_solver_module.CreateSolver(self.main_model_part, custom_settings["ale_settings"])
        print("::[ALEAdjointVMSMonolithicSolver]:: Construction finished")

    def AddVariables(self):
        super(ALEAdjointVMSMonolithicSolver, self).AddVariables()
        self.ale_solver.AddVariables()
        print("::[ALEAdjointVMSMonolithicSolver]:: Variables ADDED.")

    def AddDofs(self):
        super(ALEAdjointVMSMonolithicSolver, self).AddDofs()
        self.ale_solver.AddDofs()
        print("::[ALEAdjointVMSMonolithicSolver]:: DOFs ADDED.")

    def Initialize(self):
        super(ALEAdjointVMSMonolithicSolver, self).Initialize()
        self.ale_solver.Initialize()
        print("::[ALEAdjointVMSMonolithicSolver]:: Finished initialization.")

    def GetAdjointSolver(self):
        return super(ALEAdjointVMSMonolithicSolver, self)

    def GetMeshMotionSolver(self):
        return self.ale_solver
