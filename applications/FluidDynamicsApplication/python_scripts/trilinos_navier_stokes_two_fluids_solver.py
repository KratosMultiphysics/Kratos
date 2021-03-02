# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

# Import serial two fluid solver
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_two_fluids_solver import NavierStokesTwoFluidsSolver
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

def CreateSolver(model, custom_settings):
    return NavierStokesMPITwoFluidsSolver(model, custom_settings)

class NavierStokesMPITwoFluidsSolver(NavierStokesTwoFluidsSolver):

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)

        # Avoid using features that are not available in MPI yet
        self._bfecc_convection = self.settings["bfecc_convection"].GetBool()
        if self._bfecc_convection:
            self._bfecc_convection = False
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "BFECC is not implemented in MPI yet. Switching to standard level set convection.")

        if self.settings["formulation"].Has("surface_tension"):
            self.settings["formulation"]["surface_tension"].SetBool(False)
            self.main_model_part.ProcessInfo.SetValue(KratosFluid.SURFACE_TENSION, False)
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Surface tension is not implemented in MPI yet. Deactivating it.")

        if not self._reinitialization_type == None:
            if not self._reinitialization_type == "variational":
                self._reinitialization_type == "variational"
                KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Only variational redistancing is implemented in MPI. Switching to it.")

        if not self._distance_smoothing == None:
            if self._distance_smoothing:
                self._distance_smoothing = False
                KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Distance smoothing is not implemented in MPI yet. Deactivating it.")

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"Construction of NavierStokesMPITwoFluidsSolver finished.")

    def AddVariables(self):
        super(NavierStokesMPITwoFluidsSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"Variables for the Trilinos Two Fluid solver added correctly.")

    def ImportModelPart(self):
        ## Construct the distributed import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ImportModelPart()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"MPI model reading finished.")

    def PrepareModelPart(self):
        super(NavierStokesMPITwoFluidsSolver,self).PrepareModelPart()
        ## Construct the MPI communicators
        self.distributed_model_part_importer.CreateCommunicators()

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = KratosTrilinos.CreateCommunicator()
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
            err_msg = "Custom scheme creation is not allowed. Two-fluids Navier-Stokes elements manage the time integration internally."
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
            builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolver(
                epetra_communicator,
                guess_row_size,
                trilinos_linear_solver)
        return builder_and_solver

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosTrilinos.TrilinosNewtonRaphsonStrategy(
            computing_model_part,
            time_scheme,
            convergence_criterion,
            builder_and_solver,
            self.settings["maximum_iterations"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())

    def _CreateLevelSetConvectionProcess(self):
        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        linear_solver = self._GetLinearSolver()
        computing_model_part = self.GetComputingModelPart()
        epetra_communicator = self._GetEpetraCommunicator()

        levelset_convection_settings = self.settings["levelset_convection_settings"]
        if domain_size == 2:
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetConvectionProcess2D(
                epetra_communicator,
                computing_model_part,
                linear_solver,
                levelset_convection_settings)
        else:
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetConvectionProcess3D(
                epetra_communicator,
                computing_model_part,
                linear_solver,
                levelset_convection_settings)

        return level_set_convection_process

    def _CreateDistanceReinitializationProcess(self):
        # Construct the variational distance calculation process
        maximum_iterations = 2 #TODO: Make this user-definable
        linear_solver = self._GetLinearSolver()
        computing_model_part = self.GetComputingModelPart()
        epetra_communicator = self._GetEpetraCommunicator()
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosTrilinos.TrilinosVariationalDistanceCalculationProcess2D(
                epetra_communicator,
                computing_model_part,
                linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
        else:
            variational_distance_process = KratosTrilinos.TrilinosVariationalDistanceCalculationProcess3D(
                epetra_communicator,
                computing_model_part,
                linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        return variational_distance_process
