from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalSolver(main_model_part, custom_settings)

class StaticMechanicalSolver(BaseSolver.MechanicalSolver):
    """The solid mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    Public member variables:

    See solid_mechanics_solver.py for more information.
    """    
    def __init__(self, main_model_part, custom_settings):
        
        # Set defaults and validate custom settings.
        static_settings = KratosMultiphysics.Parameters("""
        {

        }
        """)

         # Validate and transfer settings
        self._validate_and_transfer_matching_settings(custom_settings, static_settings)
        time_integration_settings = custom_settings["time_integration_settings"]
        
        # Validate the remaining settings in the base class.
        if not time_integration_settings.Has("solution_type"): 
            time_integration_settings.AddEmptyValue("solution_type")
            time_integration_settings["solution_type"].SetString("Static") # Override defaults in the base class.
        if not time_integration_settings.Has("integration_method"): 
            time_integration_settings.AddEmptyValue("integration_method")
            time_integration_settings["integration_method"].SetString("Non-Linear") # Override defaults in the base class.

        # Construct the base solver.
        super(StaticMechanicalSolver, self).__init__(main_model_part, custom_settings)

        print("::[Static_Scheme]:: "+self.time_integration_settings["integration_method"].GetString()+" Scheme Ready")        
 
    #### Solver internal methods ####
        
    def _create_solution_scheme (self):

        integration_method   = self.time_integration_settings["integration_method"].GetString()
        
        if(integration_method == "Linear"):
            #mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            time_integration_method = KratosSolid.StaticMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementStaticScheme()            
        elif(integration_method == "Non-Linear" ):
            if(self.solving_strategy_settings["builder_type"].GetString() == "component_wise"):
                dynamic_factor = 0.0
                damp_factor_m  = 0.0
                mechanical_scheme = KratosSolid.ComponentWiseBossakScheme(damp_factor_m, dynamic_factor)
            else:
                #mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
                time_integration_method = KratosSolid.StaticMethod()
                time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
                time_integration_method.SetParameters(self.process_info)
                mechanical_scheme = KratosSolid.ResidualBasedDisplacementStaticScheme()            
 
        elif(integration_method == "RotationStatic"):
            #dynamic_factor = 0.0 
            #damp_factor_m  = 0.0
            #mechanical_scheme = KratosSolid.ResidualBasedRotationNewmarkScheme(dynamic_factor, damp_factor_m)
            time_integration_method = KratosSolid.StaticStepMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            angular_time_integration_method = KratosSolid.StaticStepRotationMethod()
            angular_time_integration_method.AddToProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, angular_time_integration_method, self.process_info)
            angular_time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementRotationStaticScheme()   
        else:
            raise Exception("Unsupported integration_method: " + integration_method)
     
        return mechanical_scheme

    def _create_mechanical_solver(self):
        if(self.solving_strategy_settings["builder_type"].GetString() == "component_wise"):
            mechanical_solver = self._create_component_wise_strategy()
        elif(self.solving_strategy_settings["line_search"].GetBool() == True):
            if(self.solving_strategy_settings["implex"].GetBool() == True):
                mechanical_solver = self._create_line_search_implex_strategy()
            else:
                mechanical_solver = self._create_line_search_strategy()
        else:
            if(self.time_integration_settings["integration_method"].GetString() == "Linear"):
                mechanical_solver = self._create_linear_strategy()
            else:
                mechanical_solver = self._create_newton_raphson_strategy()
                
        return mechanical_solver


    def _create_component_wise_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ComponentWiseNewtonRaphsonStrategy(self.model_part, 
                                                              mechanical_scheme, 
                                                              linear_solver, 
                                                              mechanical_convergence_criterion, 
                                                              builder_and_solver, 
                                                              self.solving_strategy_settings["max_iteration"].GetInt(),
                                                              self.solving_strategy_settings["compute_reactions"].GetBool(),
                                                              self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
                                                              self.solving_strategy_settings["move_mesh_flag"].GetBool())
     
    def _create_line_search_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        # KratosMultiphysics.LineSearchStrategy (alternative to test)
        return KratosSolid.ResidualBasedNewtonRaphsonLineSearchStrategy(self.model_part, 
                                                                        mechanical_scheme, 
                                                                        linear_solver, 
                                                                        mechanical_convergence_criterion, 
                                                                        builder_and_solver, 
                                                                        self.solving_strategy_settings["max_iteration"].GetInt(),
                                                                        self.solving_strategy_settings["compute_reactions"].GetBool(),
                                                                        self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
                                                                        self.solving_strategy_settings["move_mesh_flag"].GetBool())
    
    def _create_line_search_implex_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ResidualBasedNewtonRaphsonLineSearchImplexStrategy(self.model_part, 
                                                                              mechanical_scheme, 
                                                                              linear_solver, 
                                                                              mechanical_convergence_criterion, 
                                                                              builder_and_solver, 
                                                                              self.solving_strategy_settings["max_iteration"].GetInt(),
                                                                              self.solving_strategy_settings["compute_reactions"].GetBool(),
                                                                              self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
                                                                              self.solving_strategy_settings["move_mesh_flag"].GetBool())


    def _create_linear_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(self.model_part, 
                                                              mechanical_scheme, 
                                                              linear_solver, 
                                                              builder_and_solver, 
                                                              self.solving_strategy_settings["compute_reactions"].GetBool(), 
                                                              self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(), 
                                                              False, 
                                                              self.solving_strategy_settings["move_mesh_flag"].GetBool())

    
    def _create_newton_raphson_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.model_part, 
                                                                     mechanical_scheme, 
                                                                     linear_solver, 
                                                                     mechanical_convergence_criterion, 
                                                                     builder_and_solver,
                                                                     self.solving_strategy_settings["max_iteration"].GetInt(),
                                                                     self.solving_strategy_settings["compute_reactions"].GetBool(),
                                                                     self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.solving_strategy_settings["move_mesh_flag"].GetBool())

