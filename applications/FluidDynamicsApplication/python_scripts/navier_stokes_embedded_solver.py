from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
import navier_stokes_base_solver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesEmbeddedMonolithicSolver(main_model_part, custom_settings)

class NavierStokesEmbeddedMonolithicSolver(navier_stokes_base_solver.NavierStokesBaseSolver):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "embedded_solver_from_defaults",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 1.0,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "AMGCL",
                "max_iteration"       : 200,
                "tolerance"           : 1e-7,
                "provide_coordinates" : false,
                "smoother_type"       : "ilu0",
                "krylov_type"         : "lgmres",
                "coarsening_type"     : "aggregation",
                "scaling"             : true,
                "verbosity"           : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0
            },
            "periodic": "periodic",
            "move_mesh_flag": false,
            "reorder": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the element replace settings
        self._SetEmbeddedElementReplaceSettings()

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        print("Construction of NavierStokesEmbeddedSolver finished.")


    def AddVariables(self):
        ## Add base class variables
        super(NavierStokesEmbeddedMonolithicSolver, self).AddVariables()
        ## Add specific variables needed for the embedded solver
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)              # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ELEMENTAL_DISTANCES)   # Store the element nodal distance values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SOUND_VELOCITY)        # Speed of sound velocity
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)     # Nodal external pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)     # At the moment, the EmbeddedNavierStokes element works with the DYNAMIC_VISCOSITY
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_PRESSURE)          # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_VELOCITY)          # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)

        print("Monolithic embedded fluid solver variables added correctly")


    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,
                                                                            self.settings["time_order"].GetInt())

        if (self.settings["solver_type"].GetString() == "EmbeddedAusas"):
            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            self.find_nodal_neighbours_process = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part,
                                                                                               number_of_avg_elems,
                                                                                               number_of_avg_nodes)

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],   # Domain size (2,3)
                                                                                        self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1) # DOFs (3,4)

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        print ("Monolithic embedded solver initialization finished.")


    def InitializeSolutionStep(self):
        (self.bdf_process).Execute()
        if (self.settings["solver_type"].GetString() == "EmbeddedAusas"):
            (self.find_nodal_neighbours_process).Execute()
        (self.solver).InitializeSolutionStep()


    def Solve(self):
        (self.bdf_process).Execute()
        if (self.settings["solver_type"].GetString() == "EmbeddedAusas"):
            (self.find_nodal_neighbours_process).Execute()

        # Note that the first two time steps are dropped to fill the BDF buffer
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] >= 2):
            (self.solver).Solve()


    def _ExecuteAfterReading(self):
        ## Base class _ExecuteAfterReading call
        super(NavierStokesEmbeddedMonolithicSolver, self)._ExecuteAfterReading()

        ## Set the SOUND_VELOCITY value (wave velocity)
        if self.main_model_part.Properties[1].Has(KratosMultiphysics.SOUND_VELOCITY):
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = self.main_model_part.Properties[1][KratosMultiphysics.SOUND_VELOCITY]
        else:
            # If the wave velocity is not defined take a large enough value to consider the fluid as incompressible
            default_sound_velocity = 1e+12
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = default_sound_velocity
            # Set the wave velocity in the model part nodes
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.SOUND_VELOCITY, default_sound_velocity, self.main_model_part.Nodes)

        ## Set the DYNAMIC_VISCOSITY variable needed for the embedded element
        for element in self.main_model_part.Elements:
            dyn_viscosity = element.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            break

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)

        ## Construct the constitutive law needed for the embedded element
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian3DLaw()
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian2DLaw()

        ## Setting the nodal distance
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            import read_distance_from_file
            DistanceUtility = read_distance_from_file.DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            print("Distance function taken from the .mdpa input file.")
            # Recall to swap the distance sign (GiD considers d<0 in the fluid region)
            for node in self.main_model_part.Nodes:
                distance_value = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -distance_value)

    def _SetEmbeddedElementReplaceSettings(self):
        solver_type = self.settings["solver_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.settings.AddEmptyValue("element_replace_settings")

        if (solver_type == "Embedded"):
            if(domain_size == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedNavierStokes3D4N",
                    "condition_name": "NavierStokesWallCondition3D3N"
                }
                """)
            elif(domain_size == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedNavierStokes2D3N",
                    "condition_name": "NavierStokesWallCondition2D2N"
                }
                """)
            else:
                raise Exception("Domain size is not 2 or 3!!")

        elif (solver_type == "EmbeddedAusas"):
            if(domain_size == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedAusasNavierStokes3D4N",
                    "condition_name": "EmbeddedAusasNavierStokesWallCondition3D3N"
                }
                """)
            elif(domain_size == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedAusasNavierStokes2D3N",
                    "condition_name": "EmbeddedAusasNavierStokesWallCondition2D2N"
                }
                """)
            else:
                raise Exception("Domain size is not 2 or 3!!")

        elif (solver_type == "EmbeddedDevelopment"):
            if(domain_size == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedSymbolicNavierStokes3D4N",
                    "condition_name": "NavierStokesWallCondition3D3N"
                }
                """)
            elif(domain_size == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                    "element_name":"EmbeddedSymbolicNavierStokes2D3N",
                    "condition_name": "NavierStokesWallCondition2D2N"
                }
                """)
            else:
                raise Exception("Domain size is not 2 or 3!!")

        else:
            raise Exception("Wrong embedded solver type!!")
