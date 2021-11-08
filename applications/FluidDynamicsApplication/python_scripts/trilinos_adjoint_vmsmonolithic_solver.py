## Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Import applications
import KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension as KratosCFD
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication.adjoint_vmsmonolithic_solver import AdjointVMSMonolithicSolver

from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

def CreateSolver(main_model_part, custom_settings):
    return AdjointVMSMonolithicMPISolver(main_model_part, custom_settings)

class AdjointVMSMonolithicMPISolver(AdjointVMSMonolithicSolver):

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

    def AddVariables(self):
        ## Add variables from the base class
        super(self.__class__, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Variables for the AdjointVMSMonolithicMPISolver added correctly in each processor.")

    def ImportModelPart(self):
        ## Construct the distributed import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ExecutePartitioningAndReading()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "MPI model reading finished.")

    def PrepareModelPart(self):
        super(self.__class__,self).PrepareModelPart()
        ## Construct the MPI communicators
        self.distributed_model_part_importer.CreateCommunicators()

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = TrilinosApplication.CreateEpetraCommunicator(self.main_model_part.GetCommunicator().GetDataCommunicator())
        return self._epetra_communicator

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if scheme_type == "bossak":
            scheme = KratosCFD.TrilinosVelocityBossakAdjointScheme(self.settings["scheme_settings"], response_function, domain_size, domain_size + 1)
        elif scheme_type == "steady":
            scheme = KratosCFD.TrilinosSimpleSteadyAdjointScheme(response_function, domain_size, domain_size + 1)
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)
        return scheme

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return trilinos_linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateBuilderAndSolver(self):
        # Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 3:
            guess_row_size = 20*4
        else:
            guess_row_size = 10*3
        # Construct the Trilinos builder and solver
        trilinos_linear_solver = self._GetLinearSolver()
        epetra_communicator = self._GetEpetraCommunicator()
        if self.settings["consider_periodic_conditions"].GetBool():
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolverPeriodic(
                epetra_communicator,
                guess_row_size,
                trilinos_linear_solver,
                KratosFluid.PATCH_INDEX)
        else:
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(
                epetra_communicator,
                guess_row_size,
                trilinos_linear_solver)
        return builder_and_solver

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        linear_solver = self._GetLinearSolver()
        builder_and_solver = self._GetBuilderAndSolver()
        calculate_reaction_flag = False
        reform_dof_set_at_each_step = False
        calculate_norm_dx_flag = False
        move_mesh_flag = False
        return TrilinosApplication.TrilinosLinearStrategy(
            computing_model_part,
            time_scheme,
            linear_solver,
            builder_and_solver,
            calculate_reaction_flag,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag,
            move_mesh_flag)
