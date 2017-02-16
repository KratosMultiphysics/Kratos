from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application

## Checks that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import serial monolithic embedded solver
import navier_stokes_embedded_solver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesMPIEmbeddedMonolithicSolver(main_model_part, custom_settings)

class NavierStokesMPIEmbeddedMonolithicSolver(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ## Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_navier_stokes_embedded_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_GID_file",
                "distance_file_name"  : "distance_file"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": true,
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
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the element replace settings
        self.settings.AddEmptyValue("element_replace_settings")
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedNavierStokes3D4N",
                    "condition_name": "MonolithicWallCondition3D"
                }
                """)
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedNavierStokes2D3N",
                    "condition_name": "MonolithicWallCondition2D"
                }
                """)
        else:
            raise Exception("Domain size is not 2 or 3!!")

        print("Construction of NavierStokesMPIEmbeddedMonolithicSolver finished.")


    def AddVariables(self):

        super(NavierStokesMPIEmbeddedMonolithicSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMPI.mpi.world.barrier()

        if KratosMPI.mpi.rank == 0:
            print("Variables for the Trilinos monolithic embedded fluid solver added correctly.")


    def ImportModelPart(self):

        # Construct the Trilinos import model part utility
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
        # Execute the Metis partitioning and reading
        TrilinosModelPartImporter.ExecutePartitioningAndReading()
        # Call the base class execute after reading (substitute elements, set density, viscosity and constitutie law)
        super(NavierStokesMPIEmbeddedMonolithicSolver, self)._ExecuteAfterReading()
        # Call the base class set buffer size
        super(NavierStokesMPIEmbeddedMonolithicSolver, self)._SetBufferSize()
        # Construct the communicators
        TrilinosModelPartImporter.CreateCommunicators()

        if KratosMPI.mpi.rank == 0:
            print("Trilinos import model part performed correctly.")


    def AddDofs(self):

        super(NavierStokesMPIEmbeddedMonolithicSolver, self).AddDofs()
        KratosMPI.mpi.world.barrier()

        if KratosMPI.mpi.rank == 0:
            print("DOFs for the VMS Trilinos fluid solver added correctly.")


    def Initialize(self):
        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()
        
        ## If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        ## Creating the Trilinos convergence criteria
        self.conv_criteria = KratosTrilinos.TrilinosUPCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                               self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                               self.settings["relative_pressure_tolerance"].GetDouble(),
                                                               self.settings["absolute_pressure_tolerance"].GetDouble(),
                                                               self.EpetraCommunicator)

        ## Constructing the BDF process (time coefficients update)
        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,self.settings["time_order"].GetInt())

        ## Creating the Trilinos incremental update time scheme (the time integration is defined within the embedded element)
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

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
