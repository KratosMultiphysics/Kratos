from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport
import json

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    
    return FluidTransportSteadySolver(main_model_part, custom_settings)


class FluidTransportSteadySolver(object):

    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "fluid_transport_solver",
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "buffer_size": 2,  
            "echo_level":                         0,
            "clear_storage":                      false,
            "compute_reactions":                  false,
            "move_mesh_flag":                     false,
            "reform_dofs_at_each_step":           false,
            "block_builder":                      true,
            "solution_type":                      "Quasi-Static",
            "strategy_type":                      "Newton-Raphson",
            "convergence_criterion":              "And_criterion",
            "displacement_relative_tolerance":    1.0E-4,
            "displacement_absolute_tolerance":    1.0E-9,
            "residual_relative_tolerance":        1.0E-4,
            "residual_absolute_tolerance":        1.0E-9,
            "max_iteration":                      15,
            "linear_solver_settings":             {
                "solver_type":   "SuperLUSolver",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "ILU0Preconditioner",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },          
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""]
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
                
        # Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of FluidTransportSteadySolver finished")
    
    def AddVariables(self):

        ## ConvectionDiffusionSettings
        thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        thermal_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
        thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
        thermal_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
        thermal_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
        thermal_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
        thermal_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        thermal_settings.SetMeshVelocityVariable(KratosMultiphysics.MESH_VELOCITY)
        thermal_settings.SetVelocityVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, thermal_settings)

        ## Convection Variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        
        # Add thermal variables
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)



        print("Variables correctly added")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        
        for node in self.main_model_part.Nodes:
            ## Fluid dofs
            #node.AddDof(KratosMultiphysics.VELOCITY_X,KratosMultiphysics.REACTION_X)
            #node.AddDof(KratosMultiphysics.VELOCITY_Y,KratosMultiphysics.REACTION_Y)
            #node.AddDof(KratosMultiphysics.VELOCITY_Z,KratosMultiphysics.REACTION_Z)

            ## Thermal dofs
            node.AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION)
                
        print("DOFs correctly added")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            # Read ModelPart
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            
            # Create computing_model_part, set constitutive law and buffer size
            self._ExecuteAfterReading()
            
        else:
            raise Exception("Other input options are not yet implemented.")
                
        print ("Model reading finished")
    
    def Initialize(self):
        
        # Set ProcessInfo variables

        # Get the computing model parts
        self.computing_model_part = self.GetComputingModelPart()
        
        # Builder and solver creation
        builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        # Get the convergence criterion
                
        # Solver creation

        compute_norm_dx_flag = False

        self.Solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.computing_model_part, 
                                                                     scheme, 
                                                                     self.linear_solver, 
                                                                     builder_and_solver, 
                                                                     self.settings["compute_reactions"].GetBool(), 
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(), 
                                                                     compute_norm_dx_flag, 
                                                                     self.settings["move_mesh_flag"].GetBool())

        # Set echo_level
        self.Solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Solver.Check()

        print ("Initialization FluidTransportSteadySolver finished")
    
    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)
    
    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):

        if self.settings["clear_storage"].GetBool():
            self.Clear()
        
        self.Solver.Solve()


    # solve :: sequencial calls
    
    def InitializeStrategy(self):

        if self.settings["clear_storage"].GetBool():
            self.Clear()
        
        self.Solver.Initialize()

    def InitializeSolutionStep(self):
        self.Solver.InitializeSolutionStep()

    def Predict(self):
        self.Solver.Predict()

    def SolveSolutionStep(self):
        self.Solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.Solver.FinalizeSolutionStep()

    # solve :: sequencial calls

    def SetEchoLevel(self, level):
        self.Solver.SetEchoLevel(level)

    def Clear(self):
        self.Solver.Clear()
        
    def Check(self):
        self.Solver.Check()

    #### Specific internal functions ####

    def _ExecuteAfterReading(self):
        
        self.computing_model_part_name = "fluid_transport_computing_domain"
        
        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        aux_params = KratosMultiphysics.Parameters("{}")
        aux_params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)

        # CheckAndPrepareModelProcess creates the solid_computational_model_part
        import check_and_prepare_model_process_fluid_transport
        check_and_prepare_model_process_fluid_transport.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

        # # Constitutive law import
        # import poromechanics_constitutivelaw_utility
        # poromechanics_constitutivelaw_utility.SetConstitutiveLaw(self.main_model_part)

        self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        minimum_buffer_size = self.GetMinimumBufferSize()
        if(minimum_buffer_size > self.main_model_part.GetBufferSize()):
            self.main_model_part.SetBufferSize( minimum_buffer_size )

    def _ConstructBuilderAndSolver(self, block_builder):
        
        # Creating the builder and solver
        if(block_builder):
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)
        
        return builder_and_solver