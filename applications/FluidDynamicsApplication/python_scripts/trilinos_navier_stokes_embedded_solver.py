from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application

# Import serial monolithic embedded solver
import navier_stokes_embedded_solver
import trilinos_import_model_part_utility

def CreateSolver(model, custom_settings):
    return NavierStokesMPIEmbeddedMonolithicSolver(model, custom_settings)

class NavierStokesMPIEmbeddedMonolithicSolver(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver):


    def _ValidateSettings(self, settings):

        default_settings = KratosMultiphysics.Parameters("""{
            "solver_type": "Embedded",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_GID_file",
                "distance_file_name"  : "distance_file"
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"                        : "MultiLevelSolver",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-6,
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
            "move_mesh_flag": false,
            "formulation": {
                "element_type": "embedded_element_from_defaults",
                "dynamic_tau": 1.0
            },
            "fm_ale_settings": {
                "fm_ale_step_frequency": 0,
                "structure_model_part_name": "",
                "search_radius" : 1.0
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):
        # Note: deliberately calling the constructor of the base python solver (the parent of my parent)
        super(navier_stokes_embedded_solver.NavierStokesEmbeddedMonolithicSolver, self).__init__(model,custom_settings)

        self._is_printing_rank = (KratosMPI.mpi.rank == 0)

        self.min_buffer_size = 3
        self.embedded_formulation = navier_stokes_embedded_solver.EmbeddedFormulation(self.settings["formulation"])
        self.element_name = self.embedded_formulation.element_name
        self.condition_name = self.embedded_formulation.condition_name

        ## Construct the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedMonolithicSolver","Construction of NavierStokesMPIEmbeddedMonolithicSolver finished.")

    def AddVariables(self):
        super(NavierStokesMPIEmbeddedMonolithicSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedMonolithicSolver","Variables for the Trilinos monolithic embedded fluid solver added correctly.")

    def ImportModelPart(self):
        ## Construct the Trilinos import model part utility
        self.trilinos_model_part_importer = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.trilinos_model_part_importer.ImportModelPart()
        ## Sets DENSITY, VISCOSITY and SOUND_VELOCITY

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedMonolithicSolver","MPI model reading finished.")

    def PrepareModelPart(self):
        super(NavierStokesMPIEmbeddedMonolithicSolver,self).PrepareModelPart()
        ## Construct Trilinos the communicators
        self.trilinos_model_part_importer.CreateCommunicators()

    def AddDofs(self):
        super(NavierStokesMPIEmbeddedMonolithicSolver, self).AddDofs()
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedMonolithicSolver","DOFs for the VMS Trilinos fluid solver added correctly.")

    def Initialize(self):
        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

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

        # For the primitive Ausas formulation, set the find nodal neighbours process
        # Recall that the Ausas condition requires the nodal neighbouts.
        if (self.settings["formulation"]["element_type"].GetString() == "embedded_ausas_navier_stokes"):
            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            self.find_nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.GetComputingModelPart(),
                                                                                               number_of_avg_elems,
                                                                                               number_of_avg_nodes)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesMPIEmbeddedMonolithicSolver","Solver initialization finished.")