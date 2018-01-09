from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.mpi import *
KratosMultiphysics.CheckForPreviousImport()
import mesh_solver_base


def CreateSolver(model_part, custom_settings):
    return TrilinosMeshSolverBase(model_part, custom_settings)


class TrilinosMeshSolverBase(mesh_solver_base.MeshSolverBase):
    def __init__(self, model_part, custom_settings):
        super(TrilinosMeshSolverBase, self).__init__(model_part, custom_settings)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("::[TrilinosMeshSolverBase]:: Construction finished")

    #### Public user interface functions ####

    def AddVariables(self):
        super(TrilinosMeshSolverBase, self).AddVariables()
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        mpi.world.barrier()
        if mpi.rank == 0:
            print("::[MeshSolverBase]:: Variables ADDED.")

    def AddDofs(self):
        super(TrilinosMeshSolverBase, self).AddDofs()
        mpi.world.barrier()
        if mpi.rank == 0:
            print("::[MeshSolverBase]:: DOFs ADDED.")

    #### Specific internal functions ####

    def get_communicator(self):
        if not hasattr(self, '_communicator'):
            self._communicator = TrilinosApplication.CreateCommunicator()
        return self._communicator
        
    #### Private functions ####
    
    def _create_linear_solver(self):
        import trilinos_linear_elastic_ml_solver
        nit_max = 10000
        linear_tol = 1e-5
        linear_solver = trilinos_linear_elastic_ml_solver.MultilevelLinearSolver(linear_tol, nit_max)
        return linear_solver

    def _create_mesh_motion_solver(self):
        raise Exception("Mesh motion solver must be created by the derived class.")
