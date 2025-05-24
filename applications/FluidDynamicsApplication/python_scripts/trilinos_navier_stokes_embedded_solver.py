# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

# Import serial monolithic embedded solver
from KratosMultiphysics.FluidDynamicsApplication import navier_stokes_embedded_solver

from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

def CreateSolver(model, custom_settings):
    return NavierStokesMPIEmbeddedMonolithicSolver(model, custom_settings)

class NavierStokesMPIEmbeddedMonolithicSolver(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver):

    def __init__(self, model, custom_settings):
        # Call the serial base class constructor
        super().__init__(model,custom_settings)

        #TODO: Remove this once the FM-ALE is fully MPI compatible
        if self._FmAleIsActive():
            err_msg = "FM-ALE algorithm implementation is not fully MPI compatible yet. Deactivate it setting \'fm_ale_step_frequency\' to 0."
            raise Exception(err_msg)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"Construction of NavierStokesMPIEmbeddedMonolithicSolver finished.")

    def AddVariables(self):
        super(NavierStokesMPIEmbeddedMonolithicSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"Variables for the Trilinos monolithic embedded fluid solver added correctly.")

    def ImportModelPart(self):
        ## Construct the Distributed import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ImportModelPart()
        ## Sets DENSITY, VISCOSITY and SOUND_VELOCITY

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"MPI model reading finished.")

    def PrepareModelPart(self):
        super(NavierStokesMPIEmbeddedMonolithicSolver,self).PrepareModelPart()
        ## Construct MPI the communicators
        self.distributed_model_part_importer.CreateCommunicators()

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = KratosTrilinos.CreateEpetraCommunicator(self.main_model_part.GetCommunicator().GetDataCommunicator())
        return self._epetra_communicator

    def _CreateScheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Cases in which the element manages the time integration
        if self.element_integrates_in_time:
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticSchemeSlip(
                domain_size,
                domain_size + 1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        # Cases in which a time scheme manages the time integration
        else:
            err_msg = "Custom scheme creation is not allowed. Embedded Navier-Stokes elements manage the time integration internally."
            raise Exception(err_msg)
        return scheme

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return trilinos_linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateConvergenceCriterion(self):
        convergence_criterion = KratosTrilinos.TrilinosMixedGenericCriteria(
            [(KratosMultiphysics.VELOCITY, self.settings["relative_velocity_tolerance"].GetDouble(), self.settings["absolute_velocity_tolerance"].GetDouble()),
            (KratosMultiphysics.PRESSURE, self.settings["relative_pressure_tolerance"].GetDouble(), self.settings["absolute_pressure_tolerance"].GetDouble())])
        convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
        return convergence_criterion

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
            builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolverPeriodic(
                epetra_communicator,
                guess_row_size,
                trilinos_linear_solver,
                KratosFluid.PATCH_INDEX)
        else:
            if self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
                bs_params = self.settings["builder_and_solver_settings"]["advanced_settings"]
                bs_params.AddEmptyValue("guess_row_size")
                bs_params["guess_row_size"].SetInt(guess_row_size)
                builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolver(epetra_communicator, trilinos_linear_solver, bs_params)
            else:
                if self.settings["multi_point_constraints_used"].GetBool():
                    raise Exception("MPCs not yet implemented in MPI")
                if self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0:
                    KratosMultiphysics.Logger.PrintWarning("Constraints are not yet implemented in MPI and will therefore not be considered!")
                builder_and_solver = KratosTrilinos.TrilinosEliminationBuilderAndSolver(epetra_communicator, guess_row_size, trilinos_linear_solver)
        return builder_and_solver

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        linear_solver = self._GetLinearSolver()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosTrilinos.TrilinosNewtonRaphsonStrategy(
            computing_model_part,
            time_scheme,
            linear_solver,
            convergence_criterion,
            builder_and_solver,
            self.settings["maximum_iterations"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())

    @classmethod
    def _FmAleIsActive(self):
        return False