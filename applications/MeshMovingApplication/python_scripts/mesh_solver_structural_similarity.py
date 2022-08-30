# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.mesh_solver_base import MeshSolverBase


def CreateSolver(model, custom_settings):
    return MeshSolverStructuralSimilarity(model, custom_settings)


class MeshSolverStructuralSimilarity(MeshSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MeshSolverStructuralSimilarity]:: Construction finished")

    def _CreateSolutionStrategy(self):
        linear_solver = self._GetLinearSolver()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        poisson_ratio = self.settings["poisson_ratio"].GetDouble()
        solving_strategy = KratosMeshMoving.StructuralMeshMovingStrategy(self.mesh_model_part,
                                                             linear_solver,
                                                             0,
                                                             reform_dofs_each_step,
                                                             compute_reactions,
                                                             False,
                                                             self.echo_level,
                                                             poisson_ratio,
                                                             self.reinitialize_model_part_each_step)
        return solving_strategy
