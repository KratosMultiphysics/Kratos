from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver(solid_mechanics_solver.MechanicalSolver):
    
    
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
            "solver_type": "solid_mechanics_implicit_dynamic_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "implex": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]

        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of Implicit Mechanical Solver finished")
        

    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Get the solid computing model part
        self.computing_model_part = self.GetComputingModelPart()
        
        # Builder and solver creation
        builder_and_solver = self._GetBuilderAndSolver(self.settings["component_wise"].GetBool(), 
                                                       self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["scheme_type"].GetString(), 
                                                    self.settings["component_wise"].GetBool(),
                                                    self.settings["compute_contact_forces"].GetBool())
        
        # Get the convergence criterion
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        
        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     mechanical_convergence_criterion,
                                     builder_and_solver,
                                     self.settings["max_iteration"].GetInt(),
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool(),
                                     self.settings["component_wise"].GetBool(),
                                     self.settings["line_search"].GetBool(),
                                     self.settings["implex"].GetBool())

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")

        
    #### Decomposed Newton-Raphson resolution functions ####
    
    def SolverInitialize(self):
        self.mechanical_solver.Initialize()
        
    def SolverInitializeSolutionStep(self):
        self.mechanical_solver.InitializeSolutionStep()
        
    def SolverPredict(self):
        self.mechanical_solver.Predict()
        
    def SolverSolveSolutionStep(self):
        self.mechanical_solver.SolveSolutionStep()
        
    def SolverFinalizeSolutionStep(self):
        self.mechanical_solver.FinalizeSolutionStep()


    #### Specific internal functions ####
    
    def _GetSolutionScheme(self, scheme_type, component_wise, compute_contact_forces):

        if(scheme_type == "Newmark"):
            self.settings.AddEmptyValue("damp_factor_m")  
            self.settings.AddEmptyValue("dynamic_factor")
            self.settings["damp_factor_m"].SetDouble(0.0)
            self.settings["dynamic_factor"].SetDouble(1.0)
                                                                           
        elif(scheme_type == "Bossak"):
            self.settings.AddEmptyValue("damp_factor_m")  
            self.settings.AddEmptyValue("dynamic_factor")
            self.settings["damp_factor_m"].SetDouble(-0.01)
            self.settings["dynamic_factor"].SetDouble(1.0)    
        
        # Creating the implicit solution scheme:  
        if (scheme_type == "Newmark" or scheme_type == "Bossak"):
            #~ self.main_model_part.ProcessInfo[KratosSolid.RAYLEIGH_ALPHA] = 0.0
            #~ self.main_model_part.ProcessInfo[KratosSolid.RAYLEIGH_BETA ] = 0.0
          
            if(component_wise):
                mechanical_scheme = KratosSolid.ComponentWiseBossakScheme(self.settings["damp_factor_m"].GetDouble(), 
                                                                          self.settings["dynamic_factor"].GetDouble())
            else:
                mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(self.settings["damp_factor_m"].GetDouble())
                                
        return mechanical_scheme
    
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search, implex):
        if(component_wise):
            self.mechanical_solver = KratosSolid.ComponentWiseNewtonRaphsonStrategy(self.computing_model_part, 
                                                                                    mechanical_scheme, 
                                                                                    self.linear_solver, 
                                                                                    mechanical_convergence_criterion, 
                                                                                    builder_and_solver, 
                                                                                    max_iters, 
                                                                                    compute_reactions, 
                                                                                    reform_step_dofs, 
                                                                                    move_mesh_flag)
        else:
            if(line_search):
                if(implex):
                    self.mechanical_solver = KratosSolid.ResidualBasedNewtonRaphsonLineSearchImplexStrategy(self.computing_model_part, 
                                                                                                            mechanical_scheme, 
                                                                                                            self.linear_solver, 
                                                                                                            mechanical_convergence_criterion, 
                                                                                                            builder_and_solver, 
                                                                                                            max_iters, 
                                                                                                            compute_reactions, 
                                                                                                            reform_step_dofs, 
                                                                                                            move_mesh_flag)
                else:
                    self.mechanical_solver = KratosSolid.ResidualBasedNewtonRaphsonLineSearchStrategy(self.computing_model_part, 
                                                                                                      mechanical_scheme, 
                                                                                                      self.linear_solver, 
                                                                                                      mechanical_convergence_criterion, 
                                                                                                      builder_and_solver, 
                                                                                                      max_iters, 
                                                                                                      compute_reactions, 
                                                                                                      reform_step_dofs, 
                                                                                                      move_mesh_flag)


            else:
                self.mechanical_solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part, 
                                                                                               mechanical_scheme, 
                                                                                               self.linear_solver, 
                                                                                               mechanical_convergence_criterion, 
                                                                                               builder_and_solver, 
                                                                                               max_iters, 
                                                                                               compute_reactions, 
                                                                                               reform_step_dofs, 
                                                                                               move_mesh_flag)
