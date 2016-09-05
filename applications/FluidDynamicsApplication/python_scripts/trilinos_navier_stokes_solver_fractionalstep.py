from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *

# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return Trilinos_NavierStokesSolver_FractionalStep(main_model_part, custom_settings)

class Trilinos_NavierStokesSolver_FractionalStep:
    
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        if main_model_part.ProcessInfo[DOMAIN_SIZE] == 2:
            self.default_element = "FractionalStep2D3N"
            self.default_condition =  "WallCondition2D2N"
        elif main_model_part.ProcessInfo[DOMAIN_SIZE] == 3: 
            self.default_element = "FractionalStep3D4N"
            self.default_condition =  "WallCondition3D3N"
        else:
            Msg = 'Trilinos_NavierStokesSolver_FractionalStep Error:\n'
            Msg+= 'Unsupported number of dimensions: {0}\n'.format(main_model_part.ProcessInfo[DOMAIN_SIZE])
            raise Exception(Msg)
        
        ##settings string in json format
        default_settings = Parameters("""
        {
            "solver_type": "trilinos_navier_stokes_solver_fractionalstep",
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
                "solver_type"                    : "ML",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "symmetric"                      : true,
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "ML",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "symmetric"                      : false,
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[""]
        }""")
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solvers
        import trilinos_linear_solver_factory
        self.settings["pressure_linear_solver_settings"].PrettyPrintJsonString()
        self.pressure_linear_solver = trilinos_linear_solver_factory.ConstructSolver(default_settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = trilinos_linear_solver_factory.ConstructSolver(default_settings["velocity_linear_solver_settings"])

        print("Construction of Trilinos_NavierStokesSolver_FractionalStep finished")
        
    def GetMinimumBufferSize(self):
        return 3;

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(FRACT_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
        self.main_model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(CONV_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        self.main_model_part.AddNodalSolutionStepVariable(REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(PATCH_INDEX)

        mpi.world.barrier()
        
        if mpi.rank == 0:
            print("variables for the trilinos fractional step solver added correctly")
        
    def ImportModelPart(self):
        
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            
            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            self.settings["element_replace_settings"].AddEmptyValue("element_name")
            self.settings["element_replace_settings"]["element_name"].SetString(self.default_element)
            self.settings["element_replace_settings"].AddEmptyValue("condition_name")
            self.settings["element_replace_settings"]["condition_name"].SetString(self.default_condition)
                    
            ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
            
            ##here we shall check that the input read has the shape we like
            self.settings.AddEmptyValue("prepare_model_part_settings")
            self.settings["prepare_model_part_settings"].AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            self.settings["prepare_model_part_settings"].AddValue("skin_parts",self.settings["skin_parts"])

            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, self.settings["prepare_model_part_settings"]).Execute()
            
            #here we read the VISCOSITY and DENSITY and we apply it to the nodes
            for el in self.main_model_part.Elements:
                rho = el.Properties.GetValue(DENSITY)
                kin_viscosity = el.Properties.GetValue(VISCOSITY)
                break
            
            VariableUtils().SetScalarVar(DENSITY, rho, self.main_model_part.Nodes)              # Set density
            VariableUtils().SetScalarVar(VISCOSITY, kin_viscosity, self.main_model_part.Nodes)  # Set kinematic viscosity

        else:
            raise Exception("other input options are not yet implemented")
        
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("model reading finished")


    def AddDofs(self):
        
        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)
            node.AddDof(VELOCITY_X, REACTION_X)
            node.AddDof(VELOCITY_Y, REACTION_Y)
            node.AddDof(VELOCITY_Z, REACTION_Z)
        
        mpi.world.barrier()
        if mpi.rank == 0:
            print("dofs for the trilinos fractional step solver added correctly")

    
    def Initialize(self):

        self.EpetraComm = CreateCommunicator()
        
        
        compute_model_part = self.GetComputingModelPart()
        
        MoveMeshFlag = False
        
        use_slip_conditions = True

        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"] == True:
            self.solver_settings = TrilinosFractionalStepSettingsPeriodic(
                    self.EpetraComm,
                    compute_model_part,
                    compute_model_part.ProcessInfo[DOMAIN_SIZE],
                    self.settings["time_order"].GetInt(),
                    use_slip_conditions,
                    MoveMeshFlag,
                    self.settings["reform_dofs_at_each_step]"].GetBool(),
                    PATCH_INDEX
                    )
                                                                  
        else:
            self.solver_settings = TrilinosFractionalStepSettings(
                    self.EpetraComm,
                    compute_model_part,
                    compute_model_part.ProcessInfo[DOMAIN_SIZE],
                    self.settings["time_order"].GetInt(),
                    use_slip_conditions,
                    MoveMeshFlag,
                    self.settings["reform_dofs_at_each_step"].GetBool()
                    )
                                                              
        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver_settings.SetStrategy(TrilinosStrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.settings["velocity_tolerance"].GetDouble(),
                                         self.settings["maximum_velocity_iterations"].GetInt()
                                         )

        self.solver_settings.SetStrategy(TrilinosStrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.settings["pressure_tolerance"].GetDouble(),
                                         self.settings["maximum_pressure_iterations"].GetInt()
                                         )
        
        self.solver = TrilinosFSStrategy(compute_model_part,
                                         self.solver_settings,
                                         self.settings["predictor_corrector"].GetBool(),
                                         PATCH_INDEX)

        #self.solver.Check()

        print ("Initialization Trilinos_NavierStokesSolver_FractionalStep Finished")
        
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

        #if(self.compute_reactions):
        #    self.solver.CalculateReactions()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()
        
    def Check(self):
        self.solver.Check()


