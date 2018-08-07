from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication", "TrilinosApplication")

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Other imports
import KratosMultiphysics.mpi as KratosMPI
import mesh_solver_base


def CreateSolver(mesh_model_part, custom_settings):
    return TrilinosMeshSolverBase(mesh_model_part, custom_settings)


class TrilinosMeshSolverBase(mesh_solver_base.MeshSolverBase):
    def __init__(self, mesh_model_part, custom_settings):
        if not custom_settings.Has("mesh_motion_linear_solver_settings"): # Override defaults in the base class.
            linear_solver_settings = KratosMultiphysics.Parameters("""{
                "solver_type" : "AmesosSolver",
                "amesos_solver_type" : "Amesos_Klu"
            }""")
            custom_settings.AddValue("mesh_motion_linear_solver_settings", linear_solver_settings)
        super(TrilinosMeshSolverBase, self).__init__(mesh_model_part, custom_settings)
        self.print_on_rank_zero("::[TrilinosMeshSolverBase]:: Construction finished")

    #### Public user interface functions ####

    def AddVariables(self):
        super(TrilinosMeshSolverBase, self).AddVariables()
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        self.print_on_rank_zero("::[MeshSolverBase]:: Variables ADDED.")

    def AddDofs(self):
        super(TrilinosMeshSolverBase, self).AddDofs()
        self.print_on_rank_zero("::[MeshSolverBase]:: DOFs ADDED.")

    #### Specific internal functions ####

    def get_communicator(self):
        if not hasattr(self, '_communicator'):
            self._communicator = TrilinosApplication.CreateCommunicator()
        return self._communicator

    def print_on_rank_zero(self, *args):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            print(" ".join(map(str,args)))

    #### Private functions ####

    def _create_linear_solver(self):
        import trilinos_linear_solver_factory
        linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["mesh_motion_linear_solver_settings"])
        return linear_solver

    def _create_mesh_motion_solving_strategy(self):
        raise Exception("Mesh motion solver must be created by the derived class.")