from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
KratosMultiphysics.CheckForPreviousImport()
import mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return MeshSolverStructuralSimilarity(model_part, custom_settings)


class MeshSolverStructuralSimilarity(mesh_solver_base.MeshSolverBase):
    def __init__(self, model_part, custom_settings):
        super(MeshSolverStructuralSimilarity, self).__init__(model_part, custom_settings)
        print("::[MeshSolverStructuralSimilarity]:: Construction finished")

    def _create_mesh_motion_solver(self):
        linear_solver = self.get_linear_solver()
        time_order = self.settings["time_order"].GetInt()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        solver = ALEApplication.StructuralMeshMovingStrategy(self.model_part,
                                                             linear_solver,
                                                             time_order,
                                                             reform_dofs_each_step,
                                                             compute_reactions)
        return solver
