from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD     # Fluid dynamics application

# Import base class file
import navier_stokes_solver_vmsmonolithic

def CreateSolver(main_model_part, custom_settings):
    return TrilinosNavierStokesSolverMonolithic(main_model_part, custom_settings)

class TrilinosNavierStokesSolverMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic):

    def __init__(self, main_model_part, custom_settings):

        self.element_name = "VMS"
        self.condition_name = "MonolithicWallCondition"
        self.min_buffer_size = 2

        self._is_printing_rank = (KratosMPI.mpi.rank == 0)

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ## Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_navier_stokes_solver_vmsmonolithic",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
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
                "solver_type"                        : "MultiLevelSolver",
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
                "automatic_time_step" : true,
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
            "turbulence_model": "None"
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.Print("Construction of TrilinosNavierStokesSolverMonolithic finished.")


    def AddVariables(self):
        ## Add variables from the base class
        super(TrilinosNavierStokesSolverMonolithic, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.Print("Variables for the VMS fluid Trilinos solver added correctly in each processor.")


    def ImportModelPart(self):
        ## Construct the Trilinos import model part utility
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        TrilinosModelPartImporter.ExecutePartitioningAndReading()
        ## Replace default elements and conditions
        super(TrilinosNavierStokesSolverMonolithic, self)._replace_elements_and_conditions()
        ## Executes the check and prepare model process
        super(TrilinosNavierStokesSolverMonolithic, self)._execute_check_and_prepare()
        ## Call the base class set buffer size
        super(TrilinosNavierStokesSolverMonolithic, self)._set_buffer_size()
        ## Sets DENSITY, VISCOSITY and SOUND_VELOCITY
        super(TrilinosNavierStokesSolverMonolithic, self)._set_physical_properties()
        ## Construct Trilinos the communicators
        TrilinosModelPartImporter.CreateCommunicators()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.Print("MPI model reading finished.")


    def AddDofs(self):
        ## Base class DOFs addition
        super(TrilinosNavierStokesSolverMonolithic, self).AddDofs()
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.Print("DOFs for the VMS Trilinos fluid solver added correctly in all processors.")


    def Initialize(self):
        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        ## If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._get_automatic_time_stepping_utility()

        ## Creating the Trilinos convergence criteria
        self.conv_criteria = KratosTrilinos.TrilinosUPCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                               self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                               self.settings["relative_pressure_tolerance"].GetDouble(),
                                                               self.settings["absolute_pressure_tolerance"].GetDouble())

        ## Creating the Trilinos time scheme
        if (self.settings["turbulence_model"].GetString() == "None"):
            if self.settings["consider_periodic_conditions"].GetBool() == True:
                self.time_scheme = KratosTrilinos.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                          self.settings["move_mesh_strategy"].GetInt(),
                                                                                                          self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                                          KratosCFD.PATCH_INDEX)
            else:
                self.time_scheme = KratosTrilinos.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                          self.settings["move_mesh_strategy"].GetInt(),
                                                                                                          self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])


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
                                                                                           KratosCFD.PATCH_INDEX)
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

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.Print("Monolithic MPI solver initialization finished.")