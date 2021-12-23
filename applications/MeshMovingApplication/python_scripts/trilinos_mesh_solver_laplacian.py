# Importing the Kratos Library
import KratosMultiphysics

# Import applications
from KratosMultiphysics.MeshMovingApplication import TrilinosExtension as TrilinosMeshMoving

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.trilinos_mesh_solver_base import TrilinosMeshSolverBase


def CreateSolver(model, custom_settings):
    return TrilinosMeshSolverLaplacian(model, custom_settings)


class TrilinosMeshSolverLaplacian(TrilinosMeshSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        if custom_settings["buffer_size"].GetInt() < 2:
            raise Exception("A buffer_size of at least 2 is required!")
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverLaplacian]:: Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "buffer_size"           : 2
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    #### Private functions ####

    def _CreateSolutionStrategy(self):
        linear_solver = self._GetLinearSolver()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        communicator = self.get_communicator()
        solving_strategy = TrilinosMeshMoving.TrilinosLaplacianMeshMovingStrategy(
            communicator,
            self.mesh_model_part,
            linear_solver,
            0,
            reform_dofs_each_step,
            compute_reactions,
            False,
            self.echo_level)
        return solving_strategy
