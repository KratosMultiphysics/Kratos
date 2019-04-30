from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import baseclass
from KratosMultiphysics.MeshMovingApplication.trilinos_mesh_solver_base import TrilinosMeshSolverBase


def CreateSolver(model, custom_settings):
    return TrilinosMeshSolverLaplacian(model, custom_settings)


class TrilinosMeshSolverLaplacian(TrilinosMeshSolverBase):
    def __init__(self, model, custom_settings):
        if custom_settings.Has("buffer_size"):
            buffer_size = custom_settings["buffer_size"].GetInt()
            if buffer_size < 2:
                raise Exception("A buffer_size of at least 2 is required!")
        else: # overwritting baseclass-default
            custom_settings.AddEmptyValue("buffer_size").SetInt(2)
        super(TrilinosMeshSolverLaplacian, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMeshSolverLaplacian]:: Construction finished")

    #### Private functions ####

    def _create_mesh_motion_solving_strategy(self):
        linear_solver = self.get_linear_solver()
        time_order = self.settings["time_order"].GetInt()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        calculate_mesh_velocities = self.settings["calculate_mesh_velocities"].GetBool()
        echo_level = self.settings["echo_level"].GetInt()
        communicator = self.get_communicator()
        solving_strategy = TrilinosApplication.TrilinosLaplacianMeshMovingStrategy(
            communicator,
            self.mesh_model_part,
            linear_solver,
            time_order,
            reform_dofs_each_step,
            compute_reactions,
            calculate_mesh_velocities,
            echo_level)
        return solving_strategy
