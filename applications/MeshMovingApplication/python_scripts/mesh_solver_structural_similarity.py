from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Other imports
import mesh_solver_base


def CreateSolver(mesh_model_part, custom_settings):
    return MeshSolverStructuralSimilarity(mesh_model_part, custom_settings)


class MeshSolverStructuralSimilarity(mesh_solver_base.MeshSolverBase):
    def __init__(self, mesh_model_part, custom_settings):
        super(MeshSolverStructuralSimilarity, self).__init__(mesh_model_part, custom_settings)
        print("::[MeshSolverStructuralSimilarity]:: Construction finished")

    def _create_mesh_motion_solving_strategy(self):
        linear_solver = self.get_linear_solver()
        time_order = self.settings["time_order"].GetInt()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        calculate_mesh_velocities = self.settings["calculate_mesh_velocities"].GetBool()
        echo_level = self.settings["echo_level"].GetInt()
        solving_strategy = KratosMeshMoving.StructuralMeshMovingStrategy(self.mesh_model_part,
                                                             linear_solver,
                                                             time_order,
                                                             reform_dofs_each_step,
                                                             compute_reactions,
                                                             calculate_mesh_velocities,
                                                             echo_level)
        return solving_strategy
