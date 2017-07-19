from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.mpi import *
KratosMultiphysics.CheckForPreviousImport()
import trilinos_mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return TrilinosMeshSolverStructuralSimilarity(model_part, custom_settings)


class TrilinosMeshSolverStructuralSimilarity(trilinos_mesh_solver_base.TrilinosMeshSolverBase):
    def __init__(self, model_part, custom_settings):
        super(TrilinosMeshSolverStructuralSimilarity, self).__init__(model_part, custom_settings)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("::[TrilinosMeshSolverStructuralSimilarity]:: Construction finished")

    #### Private functions ####
    
    def _create_mesh_motion_solver(self):
        linear_solver = self.get_linear_solver()
        time_order = self.settings["time_order"].GetInt()
        reform_dofs_each_step = self.settings["reform_dofs_each_step"].GetBool()
        comm = self.get_communicator()
        solver = TrilinosApplication.TrilinosStructuralMeshMovingStrategy(
            comm,
            self.model_part,
            linear_solver,
            time_order,
            reform_dofs_each_step)
        return solver
