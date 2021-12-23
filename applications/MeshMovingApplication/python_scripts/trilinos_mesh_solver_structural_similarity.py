# Importing the Kratos Library
import KratosMultiphysics

# Import applications
from KratosMultiphysics.MeshMovingApplication import TrilinosExtension as TrilinosMeshMoving

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.trilinos_mesh_solver_base import TrilinosMeshSolverBase


def CreateSolver(model, custom_settings):
    return TrilinosMeshSolverStructuralSimilarity(model, custom_settings)


class TrilinosMeshSolverStructuralSimilarity(TrilinosMeshSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverStructuralSimilarity]:: Construction finished")

    #### Private functions ####

    def _CreateSolutionStrategy(self):
        linear_solver = self._GetLinearSolver()
        communicator = self.get_communicator()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        solving_strategy = TrilinosMeshMoving.TrilinosStructuralMeshMovingStrategy(
            communicator,
            self.mesh_model_part,
            linear_solver,
            0,
            reform_dofs_each_step,
            compute_reactions,
            False,
            self.echo_level)
        return solving_strategy