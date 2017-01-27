from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

#import KratosMultiphysics.MeshingApplication as KratosMeshing

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
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
            "solver_type": "navier_stokes_embedded_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
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
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

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

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        print("Construction of NavierStokesEmbeddedSolver finished.")


    def AddVariables(self):
        ## Add base class variables
        super(NavierStokesEmbeddedMonolithicSolver, self).AddVariables()
        ## Add specific variables needed for the embedded solver
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)          # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT) # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)           # Nodal volume necessary for refinement
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMeshing.ANISOTROPIC_RATIO)      # Anisotropic ratio            
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMeshing.MMG_METRIC)             # Metric for remeshing
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SOUND_VELOCITY)    # Wave velocity
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY) # At the moment, the EmbeddedNavierStokes element works with the DYNAMIC_VISCOSITY

        print("Monolithic embedded fluid solver variables added correctly")


    # def ImportModelPart(self):
    #     ## Base class ImportModelPart
    #     self._ModelPartReading()
    #     ## Replace elements and conditions, check the input reading and set KINEMATIC_VISCOSITY and DENSITY
    #     self._ExecuteAfterReading()
    #     ## Set buffer size
    #     self._SetBufferSize()
    #
    #     print ("Embedded solver model reading finished.")


    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        move_mesh_flag = False

        # Creating the solution strategy
        self.conv_criteria = KratosFluid.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                       self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                       self.settings["relative_pressure_tolerance"].GetDouble(),
                                                       self.settings["absolute_pressure_tolerance"].GetDouble())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,
                                                                            self.settings["time_order"].GetInt())

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
                                                                            move_mesh_flag)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        #~ self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        #~ self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.M, self.settings["regularization_coef"].GetDouble())

        print ("Monolithic embedded solver initialization finished.")


    def DivergenceClearance(self):
        if self.settings["divergence_clearance_steps"].GetInt() > 0:
            print("Calculating divergence-free initial condition")
            ## Initialize with a Stokes solution step
            try:
                import KratosMultiphysics.ExternalSolversApplication as KratosExternalSolvers
                smoother_type = KratosExternalSolvers.AMGCLSmoother.DAMPED_JACOBI
                solver_type = KratosExternalSolvers.AMGCLIterativeSolverType.CG
                gmres_size = 50
                max_iter = 200
                tol = 1e-7
                verbosity = 0
                stokes_linear_solver = KratosExternalSolvers.AMGCLSolver(smoother_type,
                                                                         solver_type,
                                                                         tol,
                                                                         max_iter,
                                                                         verbosity,
                                                                         gmres_size)
            except:
                pPrecond = DiagonalPreconditioner()
                stokes_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)

            stokes_process = KratosCFD.StokesInitializationProcess(self.main_model_part,
                                                                   stokes_linear_solver,
                                                                   self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                   KratosCFD.PATCH_INDEX)
            ## Copy periodic conditions to Stokes problem
            stokes_process.SetConditions(self.main_model_part.Conditions)
            ## Execute Stokes process
            stokes_process.Execute()
            stokes_process = None

            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_X, 0, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y, 0, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_Z, 0, 0.0)
    ##                vel = node.GetSolutionStepValue(VELOCITY)
    ##                for i in range(0,2):
    ##                    node.SetSolutionStepValue(VELOCITY,i,vel)

            self.settings["divergence_clearance_steps"].SetInt(0)
            print("Finished divergence clearance.")


    def SolverInitialize(self):
        self.DivergenceClearance()
        (self.solver).Initialize()


    def SolverInitializeSolutionStep(self):
        (self.bdf_process).Execute()
        (self.solver).InitializeSolutionStep()


    def Solve(self):
        self.DivergenceClearance()
        (self.bdf_process).Execute()
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
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosFluid.Newtonian3DLaw()
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosFluid.Newtonian2DLaw()

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
