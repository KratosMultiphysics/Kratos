from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as cfd
#from KratosMultiphysics.FluidDynamicsApplication import *
import KratosMultiphysics.IncompressibleFluidApplication 

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolver_VMSMonolithic(main_model_part, custom_settings)

class NavierStokesSolver_VMSMonolithic:
    
    
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
            "solver_type": "navier_stokes_solver_vmsmonolithic",
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
            "oss_switch": 0,
            "regularization_coef": 1000,
            "maximum_iterations": 30,
            "predictor_corrector": false,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "dynamic_tau": 0.0,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_iteration": true,
            "CalculateNormDxFlag": false,
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false,
            "use_des": false,
            "Cdes": 1.0,
            "wall_nodes": [],
            "spalart_allmaras_linear_solver": "None",
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "linear_solver_settings": {
                    "solver_type": "Super LU",
            }
        }""")
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solvers
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        #~ self.pressure_linear_solver = linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        #~ self.velocity_linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        print("Construction of NavierStokesSolver_VMSMonolithic finished")
        
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_TAU) # Variable stored in cfd_variables.h
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.OSS_SWITCH)  # Variable stored in cfd_variables.h
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.M)           # Variable stored in cfd_variables.h
        self.main_model_part.AddNodalSolutionStepVariable(cfd.PATCH_INDEX) # PATCH_INDEX is locally created in FluidDynamicsApp.

        # TURBULENCE MODELS ARE NOT ADDED YET
        #~ if config is not None:
            #~ if hasattr(config, "TurbulenceModel"):
                #~ if config.TurbulenceModel == "Spalart-Allmaras":
                    #~ model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
                    #~ model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
                    #~ model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
                    #~ model_part.AddNodalSolutionStepVariable(DISTANCE)
        
        print("variables for the vms fluid solver added correctly")

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
                    "element_name":"VMS3D4N",
                    "condition_name": "MonolithicWallCondition3D"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"VMS2D3N",
                    "condition_name": "Condition2D"
                    }
                    """)
            else:
                raise Exception("domain size is not 2 or 3")
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
            
            # check_skin_process is only working in 3D
            import check_and_preparemodel_process
            check_and_preparemodel_process.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()
            
            #here we read the VISCOSITY and DENSITY and we apply it to the nodes
            for el in self.main_model_part.Elements:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                nu = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                break
            
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, nu, self.main_model_part.Nodes)
            
            
            #if needed here we shall generate the constitutive laws
            #import constitutive_law_python_utility as constitutive_law_utils
            #constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(main_model_part, self.settings["DomainSize"]);
            #constitutive_law.Initialize();
        else:
            raise Exception("other input options are not yet implemented")
        
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("model reading finished")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z)

        print("dofs for the vms fluid solver added correctly")

        # TURBULENCE MODELS ARE NOT ADDED YET
        #~ if config is not None:
            #~ if hasattr(config, "TurbulenceModel"):
                #~ if config.TurbulenceModel == "Spalart-Allmaras":
                    #~ for node in model_part.Nodes:
                        #~ node.AddDof(TURBULENT_VISCOSITY)

        print("dofs for the vms monolithic solver added correctly")

    
    def Initialize(self):
        
        self.compute_model_part = self.GetComputeModelPart()
        
        MoveMeshFlag = False
        
        # PORQUE LO ASIGNABA EN EL FRACTIONAL STEP SCHEME ?¿?¿
        #~ self.settings.["use_slip_conditions"].GetBool() = True
        
        self.settings["use_slip_conditions"].SetBool(True) #ADDED TO AVOID ERROR WITH SELF.NORMAL_UTIL, IS THIS OK?
                
        # check if slip conditions are defined
        if self.settings["use_slip_conditions"].GetBool() == False:
            print('Entra1')
            for cond in self.main_model_part.Conditions:
                print('Entra2')
                print(cond.GetValue(KratosMultiphysics.IS_STRUCTURE))
                if cond.GetValue(KratosMultiphysics.IS_STRUCTURE) != 0.0:
                    print('Entra3')
                    self.settings["use_slip_conditions"] = True
                    break


                        


        # TURBULENCE MODELS ARE NOT ADDED YET
        #~ # Turbulence model
        #~ if self.use_spalart_allmaras:
            #~ self.activate_spalart_allmaras()
        
        # creating the solution strategy
        self.conv_criteria = KratosMultiphysics.IncompressibleFluidApplication.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                                                             self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                                                             self.settings["relative_pressure_tolerance"].GetDouble(),
                                                                                             self.settings["absolute_pressure_tolerance"].GetDouble())
             
                                                                                             
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.time_scheme = cfd.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                self.settings["move_mesh_strategy"].GetInt(),
                                                                                                self.compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE], 
                                                                                                cfd.PATCH_INDEX)
        else:
            self.time_scheme = cfd.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                            self.settings["move_mesh_strategy"].GetInt(),
                                                                                            self.compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
    
        # TURBULENCE MODELS ARE NOT ADDED YET
        #~ if self.turbulence_model is None:
            #~ if self.periodic == True:
                #~ self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                #~ (self.alpha, self.move_mesh_strategy, self.domain_size, PATCH_INDEX)
            #~ else:
                #~ self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                #~ (self.alpha, self.move_mesh_strategy, self.domain_size)
        #~ else:
            #~ self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                #~ (self.alpha, self.move_mesh_strategy, self.domain_size, self.turbulence_model)
            
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            builder_and_solver = cfd.ResidualBasedBlockBuilderAndSolverPeriodic(self.linear_solver, 
                                                                                cfd.PATCH_INDEX) 
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        
        
        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part, 
                                                                            self.time_scheme, 
                                                                            self.linear_solver, 
                                                                            self.conv_criteria,
                                                                            builder_and_solver, 
                                                                            self.settings["maximum_iterations"].GetInt(), 
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_iteration"].GetBool(), 
                                                                            self.settings["MoveMeshFlag"].GetBool())
                                                         
        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()
        

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.M, self.settings["regularization_coef"].GetDouble())

        print ("Initialization monolithic solver finished")
        
    def GetComputeModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString())
        
    def GetOutputVariables(self):
        pass
        
    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):

        #if self.settings["divergence_clearance_steps"].GetInt() > 0:
            #print("Calculating divergence-free initial condition")
            ## initialize with a Stokes solution step
            #try:
                #import KratosMultiphysics.ExternalSolversApplication as kes
                #smoother_type = kes.AMGCLSmoother.DAMPED_JACOBI
                #solver_type = kes.AMGCLIterativeSolverType.CG
                #gmres_size = 50
                #max_iter = 200
                #tol = 1e-7
                #verbosity = 0
                #stokes_linear_solver = kes.AMGCLSolver(
                    #smoother_type,
                    #solver_type,
                    #tol,
                    #max_iter,
                    #verbosity,
                    #gmres_size)
            #except:
                #pPrecond = DiagonalPreconditioner()
                #stokes_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)
                
            #print('HHHHHHHHHHHHOOOOOOOOOOOOLLLLLLLLLLLAAAAAAAAAAAAA')
            #err
            #stokes_process = cfd.StokesInitializationProcess(self.main_model_part,
                                                             #stokes_linear_solver,
                                                             #self.compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                             #cfd.PATCH_INDEX)
            ## copy periodic conditions to Stokes problem
            #stokes_process.SetConditions(self.main_model_part.Conditions)
            ## execute Stokes process
            #stokes_process.Execute()
            #stokes_process = None

            #for node in self.main_model_part.Nodes:
                #node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                #node.SetSolutionStepValue(ACCELERATION_X, 0, 0.0)
                #node.SetSolutionStepValue(ACCELERATION_Y, 0, 0.0)
                #node.SetSolutionStepValue(ACCELERATION_Z, 0, 0.0)
##                vel = node.GetSolutionStepValue(VELOCITY)
##                for i in range(0,2):
##                    node.SetSolutionStepValue(VELOCITY,i,vel)

            #self.settings["divergence_clearance_steps"] = 0
            #print("Finished divergence clearance")


        #if self.settings["reform_dofs_at_each_iteration"]:
            #if self.settings["use_slip_conditions"]:
                #self.normal_util.CalculateOnSimplex(self.main_model_part,
                                                    #self.compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                    #KratosMultiphysics.IS_STRUCTURE)
            #~ if self.settings["use_spalart_allmaras"]:
                #~ KratosMultiphysics.neighbour_search.Execute()
        
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()
        
    def Check(self):
        self.solver.Check()


