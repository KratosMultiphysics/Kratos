from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD     # Fluid dynamics application
from KratosMultiphysics.FluidDynamicsApplication import TrilinosExtension as TrilinosFluid
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication import navier_stokes_solver_vmsmonolithic
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility
from KratosMultiphysics.FluidDynamicsApplication.turbulence_model_solver import CreateTurbulenceModel

def CreateSolver(model, custom_settings):
    return TrilinosNavierStokesSolverMonolithic(model, custom_settings)

class TrilinosNavierStokesSolverMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic):

    @classmethod
    def GetDefaultSettings(cls):
        ## Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_navier_stokes_solver_vmsmonolithic",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "vms"
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
            "linear_solver_settings"       : {
                "solver_type"                        : "multi_level",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-8,
                "max_levels"                         : 3,
                "symmetric"                          : false,
                "reform_preconditioner_at_each_step" : true,
                "scaling"                            : true
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "time_scheme": "bossak",
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "move_mesh_flag": false,
            "turbulence_model_solver_settings": {}
        }""")

        default_settings.AddMissingParameters(super(TrilinosNavierStokesSolverMonolithic, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        # Note: deliberately calling the constructor of the base python solver (the parent of my parent)
        custom_settings = self._BackwardsCompatibilityHelper(custom_settings)
        super(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic, self).__init__(model,custom_settings)

        self.formulation = navier_stokes_solver_vmsmonolithic.StabilizedFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties

        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bossak":
            self.min_buffer_size = 2
        elif scheme_type == "bdf2":
            self.min_buffer_size = 3
        elif scheme_type == "steady":
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"bossak\", \"bdf2\" or \"steady\".\n"
            raise Exception(msg)

        ## Construct the turbulence model solver
        if not self.settings["turbulence_model_solver_settings"].IsEquivalentTo(KratosMultiphysics.Parameters("{}")):
            self._turbulence_model_solver = CreateTurbulenceModel(model, self.settings["turbulence_model_solver_settings"])
            self.condition_name = self._turbulence_model_solver.GetFluidVelocityPressureConditionName()
            KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverMonolithic", "Using " + self.condition_name + " as wall condition")

        KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverMonolithic", "Construction of TrilinosNavierStokesSolverMonolithic finished.")

    def AddVariables(self):
        ## Add variables from the base class
        super(TrilinosNavierStokesSolverMonolithic, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMultiphysics.Logger.Print("Variables for the VMS fluid Trilinos solver added correctly in each processor.")

    def ImportModelPart(self):
        ## Construct the MPI import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ImportModelPart()

        KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverMonolithic","MPI model reading finished.")

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self._SetPhysicalProperties()

        super(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic, self).PrepareModelPart()

        self.distributed_model_part_importer.CreateCommunicators()

        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.PrepareModelPart()

    def AddDofs(self):
        ## Base class DOFs addition
        super(TrilinosNavierStokesSolverMonolithic, self).AddDofs()

        KratosMultiphysics.Logger.Print("DOFs for the VMS Trilinos fluid solver added correctly in all processors.")

    def Initialize(self):
        # If there is turbulence modelling, set the Epetra communicator in the turbulence solver
        # TODO: Check if this needs to be done before creating the strategy
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.SetCommunicator(self.get_epetra_communicator())

        # Construct and set the solution strategy
        solution_strategy = self.get_solution_strategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize the solution strategy
        if not self.is_restarted():
            # Set the formulation instance process info variables
            self.formulation.SetProcessInfo(self.GetComputingModelPart())
            # Initialize the solution strategy
            solution_strategy.Initialize()
        else:
            # This try is required in case SetInitializePerformedFlag is not a member of the strategy
            try:
                solution_strategy.SetInitializePerformedFlag(True)
            except AttributeError:
                pass

        # If there is turbulence modelling, set the new solution strategy as parent strategy
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.SetParentSolvingStrategy(solution_strategy)

        KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverMonolithic", "Solver initialization finished.")

    def Finalize(self):
        self.get_solution_strategy().Clear()

        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.Finalize()

    def get_epetra_communicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = KratosTrilinos.CreateCommunicator()
        return self._epetra_communicator

    def _create_solution_scheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Cases in which the element manages the time integration
        if self.element_integrates_in_time:
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            solution_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticSchemeSlip(
                domain_size,
                domain_size + 1)
            # In case the BDF2 scheme is used inside the element, set the time discretization utility to compute the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = self.settings["time_order"].GetInt()
                if time_order == 2:
                    self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
                else:
                    raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))
            else:
                err_msg = "Requested elemental time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        # Cases in which a time scheme manages the time integration
        else:
            if not hasattr(self, "_turbulence_model_solver"):
                if self.settings["time_scheme"].GetString() == "bossak":
                    # TODO: Can we remove this periodic check, Is the PATCH_INDEX used in this scheme?
                    if self.settings["consider_periodic_conditions"].GetBool() == True:
                        solution_scheme = TrilinosFluid.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            domain_size,
                            KratosCFD.PATCH_INDEX)
                    else:
                        solution_scheme = TrilinosFluid.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.settings["move_mesh_strategy"].GetInt(),
                            domain_size)
                elif self.settings["time_scheme"].GetString() == "steady":
                    solution_scheme = TrilinosFluid.TrilinosResidualBasedSimpleSteadyScheme(
                            self.settings["velocity_relaxation"].GetDouble(),
                            self.settings["pressure_relaxation"].GetDouble(),
                            domain_size)
            else:
                self._turbulence_model_solver.Initialize()
                if self.settings["time_scheme"].GetString() == "bossak":
                    solution_scheme = TrilinosFluid.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.settings["move_mesh_strategy"].GetInt(),
                            domain_size,
                            self.settings["turbulence_model_solver_settings"]["velocity_pressure_relaxation_factor"].GetDouble(),
                            self._turbulence_model_solver.GetTurbulenceSolvingProcess())
                # Time scheme for steady state fluid solver
                elif self.settings["time_scheme"].GetString() == "steady":
                    solution_scheme = TrilinosFluid.TrilinosResidualBasedSimpleSteadyScheme(
                            self.settings["velocity_relaxation"].GetDouble(),
                            self.settings["pressure_relaxation"].GetDouble(),
                            domain_size,
                            self._turbulence_model_solver.GetTurbulenceSolvingProcess())
        return solution_scheme

    def _create_linear_solver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return trilinos_linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _create_convergence_criterion(self):
        if self.settings["time_scheme"].GetString() == "steady":
            convergence_criterion = KratosTrilinos.TrilinosResidualCriteria(
                self.settings["relative_velocity_tolerance"].GetDouble(),
                self.settings["absolute_velocity_tolerance"].GetDouble())
        else:
            convergence_criterion = KratosTrilinos.TrilinosUPCriteria(
                self.settings["relative_velocity_tolerance"].GetDouble(),
                self.settings["absolute_velocity_tolerance"].GetDouble(),
                self.settings["relative_pressure_tolerance"].GetDouble(),
                self.settings["absolute_pressure_tolerance"].GetDouble())
        convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
        return convergence_criterion

    def _create_builder_and_solver(self):
        # Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 3:
            guess_row_size = 20*4
        else:
            guess_row_size = 10*3
        # Construct the Trilinos builder and solver
        trilinos_linear_solver = self.get_linear_solver()
        epetra_communicator = self.get_epetra_communicator()
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

    def _create_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
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