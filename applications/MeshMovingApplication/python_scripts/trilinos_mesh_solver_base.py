from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication", "TrilinosApplication")

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
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
        self.print_on_rank_zero("::[TrilinosMeshSolverBase]:: Variables ADDED.")

    def ImportModelPart(self):
        self.print_on_rank_zero("::[TrilinosMeshSolverBase]:: ", "Importing model part.")
        from trilinos_import_model_part_utility import TrilinosImportModelPartUtility
        self.trilinos_model_part_importer = TrilinosImportModelPartUtility(self.mesh_model_part, self.settings)
        self.trilinos_model_part_importer.ImportModelPart()
        self.print_on_rank_zero("::[TrilinosMeshSolverBase]:: ", "Finished importing model part.")

    def PrepareModelPart(self):
        super(TrilinosMeshSolverBase, self).PrepareModelPart()
        # Construct the mpi-communicator
        self.trilinos_model_part_importer.CreateCommunicators()
        self.print_on_rank_zero("::[TrilinosMeshSolverBase]::", "ModelPart prepared for Solver.")

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