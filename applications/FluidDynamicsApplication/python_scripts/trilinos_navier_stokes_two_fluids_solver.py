from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application

# Import serial two fluid solver
import navier_stokes_two_fluids_solver
import trilinos_import_model_part_utility

def CreateSolver(model, custom_settings):
    return NavierStokesMPITwoFluidsSolver(model, custom_settings)

class NavierStokesMPITwoFluidsSolver(navier_stokes_two_fluids_solver.NavierStokesTwoFluidsSolver):

    def _ValidateSettings(self, settings):

        default_settings = KratosMultiphysics.Parameters("""{
            "solver_type": "TwoFluids",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "consider_periodic_conditions"	: false,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings":   {
                "solver_type": "AmgclMPISolver",
                "tolerance": 1.0e-6,
                "max_iteration": 10,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "None",
                "krylov_type": "cg"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "move_mesh_flag": false,
            "formulation": {
                "dynamic_tau": 1.0
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):
        # the constructor of the "grand-parent" (jumping constructor of parent) is called to avoid conflicts in attribute settings
        super(navier_stokes_two_fluids_solver.NavierStokesTwoFluidsSolver, self).__init__(model,custom_settings)

        self.element_name = "TwoFluidNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.min_buffer_size = 3

        self._is_printing_rank = (KratosMPI.mpi.rank == 0)

        if (self.settings["solver_type"].GetString() == "TwoFluids"):
            self.element_name = "TwoFluidNavierStokes"

        ## Construct the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPITwoFluidsSolver","Construction of NavierStokesMPITwoFluidsSolver finished.")


    def AddVariables(self):

        super(NavierStokesMPITwoFluidsSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPITwoFluidsSolver","Variables for the Trilinos Two Fluid solver added correctly.")

    def ImportModelPart(self):
        ## Construct the Trilinos import model part utility
        self.trilinos_model_part_importer = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.trilinos_model_part_importer.ImportModelPart()
        ## Sets DENSITY, VISCOSITY and SOUND_VELOCITY

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPITwoFluidsSolver","MPI model reading finished.")

    def PrepareModelPart(self):
        super(NavierStokesMPITwoFluidsSolver,self).PrepareModelPart()
        ## Construct Trilinos the communicators
        self.trilinos_model_part_importer.CreateCommunicators()

    def AddDofs(self):
        super(NavierStokesMPITwoFluidsSolver, self).AddDofs()
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPITwoFluidsSolver","DOFs for the Trilinos Two Fluid solver added correctly.")


    def Initialize(self):

        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.computing_model_part, self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part, 10, 10)
        (self.neighbour_search).Execute()

        self.accelerationLimitationUtility = KratosMultiphysics.FluidDynamicsApplication.AccelerationLimitationUtilities(self.computing_model_part, 5.0)

        ## If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        ## Creating the Trilinos convergence criteria
        self.conv_criteria = KratosTrilinos.TrilinosUPCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                               self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                               self.settings["relative_pressure_tolerance"].GetDouble(),
                                                               self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        ## Constructing the BDF process (time coefficients update)
        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,self.settings["time_order"].GetInt())

        #### ADDING NEW PROCESSES : level-set-convection and variational-distance-process
        self.level_set_convection_process = self._set_level_set_convection_process()
        self.variational_distance_process = self._set_variational_distance_process()

        ## Creating the Trilinos incremental update time scheme (the time integration is defined within the TwoFluidNavierStokes element)
        self.time_scheme = KratosTrilinos.TrilinosResidualBasedIncrementalUpdateStaticSchemeSlip(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],   # Domain size (2,3)
                                                                                                 self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1) # DOFs (3,4)


        ## Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            guess_row_size = 20*4
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            guess_row_size = 10*3

        ## Construct the Trilinos builder and solver
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolverPeriodic(self.EpetraCommunicator,
                                                                                           guess_row_size,
                                                                                           self.trilinos_linear_solver,
                                                                                           KratosFluid.PATCH_INDEX)
        else:
            self.builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                   guess_row_size,
                                                                                   self.trilinos_linear_solver)

        ## Construct the Trilinos Newton-Raphson strategy
        self.solver = KratosTrilinos.TrilinosNewtonRaphsonStrategy(self.main_model_part,
                                                                   self.time_scheme,
                                                                   self.trilinos_linear_solver,
                                                                   self.conv_criteria,
                                                                   self.builder_and_solver,
                                                                   self.settings["maximum_iterations"].GetInt(),
                                                                   self.settings["compute_reactions"].GetBool(),
                                                                   self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                   self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        (self.solver).Initialize()
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["formulation"]["dynamic_tau"].GetDouble())


    def _set_level_set_convection_process(self):
        # Construct the level set convection process
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetConvectionProcess2D(
                self.EpetraCommunicator,
                KratosMultiphysics.DISTANCE,
                self.computing_model_part,
                self.trilinos_linear_solver)
        else:
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetConvectionProcess3D(
                self.EpetraCommunicator,
                KratosMultiphysics.DISTANCE,
                self.computing_model_part,
                self.trilinos_linear_solver)

        return level_set_convection_process


    def _set_variational_distance_process(self):
        # Construct the variational distance calculation process
        maximum_iterations = 2
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosTrilinos.TrilinosVariationalDistanceCalculationProcess2D(
                self.EpetraCommunicator,
                self.computing_model_part,
                self.trilinos_linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
        else:
            variational_distance_process = KratosTrilinos.TrilinosVariationalDistanceCalculationProcess3D(
                self.EpetraCommunicator,
                self.computing_model_part,
                self.trilinos_linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        return variational_distance_process