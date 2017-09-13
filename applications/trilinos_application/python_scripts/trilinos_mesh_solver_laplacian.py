from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.mpi import *
KratosMultiphysics.CheckForPreviousImport()
import trilinos_mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return TrilinosMeshSolverComponentwise(model_part, custom_settings)


class TrilinosMeshSolverComponentwise(trilinos_mesh_solver_base.TrilinosMeshSolverBase):
    def __init__(self, model_part, custom_settings):
        super(TrilinosMeshSolverComponentwise, self).__init__(model_part, custom_settings)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("::[TrilinosMeshSolverComponentwise]:: Construction finished")

    #### Private functions ####

    def _create_linear_solver(self):
        import PressureMultiLevelSolver
        pressure_nit_max = 1000
        pressure_linear_tol = 1e-6
        linear_solver = PressureMultiLevelSolver.MultilevelLinearSolver(pressure_linear_tol, pressure_nit_max)
        return linear_solver

    def _create_mesh_motion_solver(self):
        domain_size = self.model_part.ProcessInfo[DOMAIN_SIZE]
        linear_solver = self.get_linear_solver()
        time_order = self.settings("time_order").GetInt()
        reform_dofs_each_step = self.settings("reform_dofs_each_step").GetBool()
        comm = self.get_communicator()
        solver = TrilinosApplication.TrilinosLaplacianMeshMovingStrategy(
            comm,
            self.model_part, 
            linear_solver, 
            domain_size,
            time_order, 
            reform_dof_at_every_step)
        return solver
