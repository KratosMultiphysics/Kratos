from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolver_FractionalStep(main_model_part, custom_settings)

class NavierStokesSolver_FractionalStep:
    
    
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
            "solver_type": "navier_stokes_solver_fractionalstep",
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "predictor_corrector": false,
            "maximum_velocity_iterations": 3,
            "maximum_pressure_iterations": 3,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "echo_level": 1,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "dynamic_tau": 0.001,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": false,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "AMGCL",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "damped_jacobi",
                "krylov_type"                    : "cg",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "AMGCL",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "bicgstab",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[""],
            "no_skin_parts":[""]
        }""")
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solvers
        import linear_solver_factory
        self.pressure_linear_solver = linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        print("Construction of NavierStokesSolver_FractionalStep finished.")
        
    def GetMinimumBufferSize(self):
        return 3;

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FRACT_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE_OLD_IT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESS_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONV_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE) #TODO: verify if it is actually used
        # Stokes needs it (in case periodic conditions are required)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)
        
        print("Variables for the VMS fluid solver added correctly.")

    def ImportModelPart(self):
        
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            
            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"FractionalStep3D4N",
                    "condition_name": "WallCondition3D3N"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"FractionalStep2D3N",
                    "condition_name": "WallCondition2D2N"
                    }
                    """)
            else:
                raise Exception("Domain size is not 2 or 3.")
        
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
            
            ## Auxiliar parameters to execute the CheckAndPrepareModelProcess
            prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
            prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])
            
            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()
            
            ## Read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            for el in self.main_model_part.Elements:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                break
            
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
            
        else:
            raise Exception("Other input options are not yet implemented.")
        
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("Model reading finished.")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z)

        print("DOFs for the VMS fluid solver added correctly")

    
    def Initialize(self):
        
        
        compute_model_part = self.GetComputingModelPart()
        
        MoveMeshFlag = False
        
        self.use_slip_conditions = True

        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"] == True:
            self.solver_settings = KratosCFD.FractionalStepSettingsPeriodic(compute_model_part,
                                                                            compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                            self.settings.GetInt(),
                                                                            self.use_slip_conditions,
                                                                            MoveMeshFlag,
                                                                            self.settings["reform_dofs_at_each_step]"].GetBool(),
                                                                            KratosCFD.PATCH_INDEX)
                                                                  
        else:
            self.solver_settings = KratosCFD.FractionalStepSettings(compute_model_part,
                                                                    compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                    self.settings["time_order"].GetInt(),
                                                                    self.use_slip_conditions,
                                                                    MoveMeshFlag,
                                                                    self.settings["reform_dofs_at_each_step"].GetBool())
                                                              
        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver_settings.SetStrategy(KratosCFD.StrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.settings["velocity_tolerance"].GetDouble(),
                                         self.settings["maximum_velocity_iterations"].GetInt())

        self.solver_settings.SetStrategy(KratosCFD.StrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.settings["pressure_tolerance"].GetDouble(),
                                         self.settings["maximum_pressure_iterations"].GetInt())


        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.solver = KratosCFD.FSStrategy(compute_model_part,
                                               self.solver_settings, 
                                               self.settings["predictor_corrector"].GetBool(),
                                               KratosCFD.PATCH_INDEX)
        else:
            self.solver = KratosCFD.FSStrategy(compute_model_part,
                                               self.solver_settings, 
                                               self.settings["predictor_corrector"].GetBool())

        self.solver.Check()
        
        print ("Initialization NavierStokesSolver_FractionalStep finished.")
        
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
    
    def Solve(self):
        self.solver.Solve()

        if(self.compute_reactions):
            self.solver.CalculateReactions()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()
        
    def Check(self):
        self.solver.Check()
