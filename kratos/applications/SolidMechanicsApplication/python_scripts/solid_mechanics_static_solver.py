from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import solid_mechanics_implicit_dynamic_solver

def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalSolver(main_model_part, custom_settings)

class StaticMechanicalSolver(solid_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    
    
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
            "solver_type": "solid_mechanics_static_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "echo_level": 0,
            "time_integration_method": "Implicit",
            "analysis_type": "nonlinear",
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_iteration": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "solution_type": "Static",
            "scheme_type": "Bossak",
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-4,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "Super LU",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "processes_sub_model_part_list": [""],
            "problem_domain_sub_model_part": "solid_model_part"
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of MechanicalSolver finished")
    
    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Builder and solver creation
        builder_and_solver = self._GetBuilderAndSolver(self.settings["component_wise"].GetBool(), 
                                                       self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["analysis_type"].GetString(), 
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
                                     self.settings["reform_dofs_at_each_iteration"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool(),
                                     self.settings["component_wise"].GetBool(),
                                     self.settings["line_search"].GetBool())

        # Set the stabilization factor
        self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")
        
    #### Specific internal functions ####
        
    def _GetSolutionScheme(self, analysis_type, component_wise, compute_contact_forces):

        if(analysis_type == "linear"):
            mechanical_scheme = SolidMechanicsApplication.ResidualBasedStaticScheme()
            
        if(analysis_type == "nonlinear" ):
            self.settings.AddEmptyValue("damp_factor_m")  
            self.settings.AddEmptyValue("dynamic_factor")
            self.settings["damp_factor_m"].SetDouble(0.0)
            self.settings["dynamic_factor"].SetDouble(0.0) # Quasi-static scheme
            
            if component_wise:
                mechanical_scheme = SolidMechanicsApplication.ComponentWiseBossakScheme(self.settings["damp_factor_m"].GetDouble(), 
                                                                                        self.settings["dynamic_factor"].GetDouble())
            else:
                if compute_contact_forces:
                    mechanical_scheme = SolidMechanicsApplication.ResidualBasedContactBossakScheme(self.settings["damp_factor_m"].GetDouble(), 
                                                                                                   self.settings["dynamic_factor"].GetDouble())
                else:
                    mechanical_scheme = SolidMechanicsApplication.ResidualBasedBossakScheme(self.settings["damp_factor_m"].GetDouble(), 
                                                                                            self.settings["dynamic_factor"].GetDouble())
                                
        return mechanical_scheme
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search):
        
        if(component_wise):
            self.mechanical_solver = SolidMechanicsApplication.ComponentWiseNewtonRaphsonStrategy(self.main_model_part, 
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
                self.mechanical_solver = SolidMechanicsApplication.ResidualBasedNewtonRaphsonLineSearchStrategy(self.main_model_part, 
                                                                                                                mechanical_scheme, 
                                                                                                                self.linear_solver, 
                                                                                                                mechanical_convergence_criterion, 
                                                                                                                builder_and_solver, 
                                                                                                                max_iters, 
                                                                                                                compute_reactions, 
                                                                                                                reform_step_dofs, 
                                                                                                                move_mesh_flag)

            else:
                if self.settings["analysis_type"].GetString() == "linear":
                    self.mechanical_solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part, 
                                                                                            mechanical_scheme, 
                                                                                            self.linear_solver, 
                                                                                            builder_and_solver, 
                                                                                            compute_reactions, 
                                                                                            reform_step_dofs, 
                                                                                            False, 
                                                                                            move_mesh_flag)
                
                else:
                    self.mechanical_solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part, 
                                                                                                   mechanical_scheme, 
                                                                                                   self.linear_solver, 
                                                                                                   mechanical_convergence_criterion, 
                                                                                                   builder_and_solver, 
                                                                                                   max_iters, 
                                                                                                   compute_reactions, 
                                                                                                   reform_step_dofs, 
                                                                                                   move_mesh_flag)
