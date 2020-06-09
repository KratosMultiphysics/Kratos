from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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
        super(MeshSolverStructuralSimilarity, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MeshSolverStructuralSimilarity]:: Construction finished")

    def _create_mesh_motion_solving_strategy(self):
        linear_solver = self.get_linear_solver()
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
                                                             poisson_ratio)
        return solving_strategy
