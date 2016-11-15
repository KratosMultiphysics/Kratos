from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesEmbeddedMonolithicSolver(main_model_part, custom_settings)

class NavierStokesEmbeddedMonolithicSolver:


    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part. This is needed since at the point of constructing the
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the model is already filled
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
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "oss_switch": 0,
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
                "solver_type"         : "AMGCL",
                "max_iteration"       : 200,
                "tolerance"           : 1e-8,
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
            "alpha":-0.1,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false
        }""")

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of NavierStokesEmbeddedSolver finished.")

    def GetMinimumBufferSize(self):
        return 3;

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_TAU) 
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.OSS_SWITCH)  
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.M)           
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)           
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.PATCH_INDEX)

        print("Variables for the fluid embedded solver added correctly.")


    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            ##here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])

            #here we replace the dummy elements we read with proper elements
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

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()
            
            #here we read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            for el in self.main_model_part.Elements:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
                break

            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)
                        
            # Constitutive law import
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosFluid.Newtonian3DLaw()
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosFluid.Newtonian2DLaw()

        else:
            raise Exception("Other input options are not yet implemented.")

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("Model reading finished.")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # Adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z)

        print("DOFs for the embedded fluid solver added correctly.")


    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        move_mesh_flag = False

        # Creating the solution strategy
        self.conv_criteria = KratosFluid.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                       self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                       self.settings["relative_pressure_tolerance"].GetDouble(),
                                                       self.settings["absolute_pressure_tolerance"].GetDouble())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,2)
        
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
        
        self.solver.Initialize()
        self.solver.Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        #~ self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        #~ self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.M, self.settings["regularization_coef"].GetDouble())

        print ("Monolithic embedded solver initialization finished.")

    def GetComputingModelPart(self):
        # Get as computational model part the "volume_model_part_name" in the ProjectParameters Json string
        #~ return self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString())

        # Get as computational model part the submodelpart generated in CheckAndPrepareModelProcess
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def DivergenceClearance(self):

        if self.settings["divergence_clearance_steps"].GetInt() > 0:
            print("Calculating divergence-free initial condition...")
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

            stokes_process = KratosFluid.StokesInitializationProcess(self.main_model_part,
                                                                     stokes_linear_solver,
                                                                     self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                     KratosFluid.PATCH_INDEX)
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
        self.solver.Initialize()

    def SolverInitializeSolutionStep(self):
        self.bdf_process.Execute()
        self.solver.InitializeSolutionStep()

    def SolverPredict(self):
        self.solver.Predict()

    def SolverSolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def Solve(self):
        self.bdf_process.Execute()
        self.DivergenceClearance()
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()
