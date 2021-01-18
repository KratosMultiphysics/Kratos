# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.mesh_solver_base import MeshSolverBase


def CreateSolver(model, custom_settings):
    return MeshSolverLaplacian(model, custom_settings)


class MeshSolverLaplacian(MeshSolverBase):
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        if custom_settings["buffer_size"].GetInt() < 2:
            raise Exception("A buffer_size of at least 2 is required!")
        KM.Logger.PrintInfo("::[MeshSolverLaplacian]:: Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "buffer_size"           : 2
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def _create_mesh_motion_solving_strategy(self):
        linear_solver = self.get_linear_solver()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        solving_strategy = KratosMeshMoving.LaplacianMeshMovingStrategy(self.mesh_model_part,
                                                            linear_solver,
                                                            0,
                                                            reform_dofs_each_step,
                                                            compute_reactions,
                                                            False,
                                                            self.echo_level)
        return solving_strategy



