from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver(BaseSolver.MechanicalSolver):
    """The solid mechanics implicit dynamic solver.

    This class creates the mechanical solvers for implicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See solid_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings): 
        
        # Set defaults and validate custom settings.
        ##TODO : solving_strategy_settings must be time_integration_settings (GiD interface changes needed)
        implicit_solver_settings = KratosMultiphysics.Parameters("""
        {
            "solving_strategy_settings":{ 
                "bossak_factor" :-0.3,
                "dynamic_factor": 1.0,
                "lumped_mass_matrix" : true,
                "consistent_mass_matrix" : false,
                "rayleigh_damping": false, 
                "rayleigh_alpha": 0.0,
                "rayleigh_beta" : 0.0
            }
        }
        """)

        # Validate and transfer settings
        if( custom_settings.Has("solving_strategy_settings") ):
            self._validate_and_transfer_matching_settings(custom_settings["solving_strategy_settings"], implicit_solver_settings["solving_strategy_settings"])
        self.implicit_solver_settings = implicit_solver_settings["solving_strategy_settings"]
        
        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)
        
        print("::[Implicit_Scheme]:: "+self.time_integration_settings["integration_method"].GetString()+" Scheme Ready")


    def GetVariables(self):

        nodal_variables = super(ImplicitMechanicalSolver, self).GetVariables()
        if(self.solving_strategy_settings["builder_type"].GetString() == "component_wise"):
            nodal_variables = nodal_variables + ['INTERNAL_FORCE','EXTERNAL_FORCE']
                
        return nodal_variables
        
    #### Solver internal methods ####
    
    def _create_solution_scheme(self):

        integration_method   = self.time_integration_settings["integration_method"].GetString()
        
        if( self.implicit_solver_settings["rayleigh_damping"].GetBool() == True ):
            self.process_info[KratosSolid.RAYLEIGH_ALPHA] = self.implicit_solver_settings["rayleigh_alpha"].GetDouble()
            self.process_info[KratosSolid.RAYLEIGH_BETA]  = self.implicit_solver_settings["rayleigh_beta"].GetDouble()
        else:
            self.process_info[KratosSolid.RAYLEIGH_ALPHA] = 0.0
            self.process_info[KratosSolid.RAYLEIGH_BETA]  = 0.0

        # compute mass lumped matrix
        if( self.implicit_solver_settings["lumped_mass_matrix"].GetBool() == True ):
            self.process_info[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = True
        else:
            # compute consistent dynamic tangent/mass matrix
            if( self.implicit_solver_settings["consistent_mass_matrix"].GetBool() == True ):
                self.process_info[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT] = True 
 
        if(self.solving_strategy_settings["builder_type"].GetString() == "component_wise"):
            dynamic_factor = self.implicit_solver_settings["dynamic_factor"].GetDouble()        
            damp_factor_m  = self.implicit_solver_settings["bossak_factor"].GetDouble()
            mechanical_scheme = KratosSolid.ComponentWiseBossakScheme(damp_factor_m, dynamic_factor)
        elif(integration_method == "Newmark"):
            #damp_factor_m = 0.0
            #mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
            time_integration_method = KratosSolid.NewmarkMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementNewmarkScheme()            
        elif(integration_method == "Bossak"):
            bossak_factor = self.implicit_solver_settings["bossak_factor"].GetDouble()
            self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = bossak_factor;
            #mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(bossak_factor)
            time_integration_method = KratosSolid.BossakMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementBossakScheme()
        elif(integration_method == "Simo"):
            bossak_factor = self.implicit_solver_settings["bossak_factor"].GetDouble()
            self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = bossak_factor;
            #mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(bossak_factor)
            time_integration_method = KratosSolid.SimoMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementSimoScheme()            
        elif(integration_method == "RotationNewmark"):
            #dynamic_factor = self.implicit_solver_settings["dynamic_factor"].GetDouble() # 0,1 
            #damp_factor_m = self.implicit_solver_settings["bossak_factor"].GetDouble()
            #mechanical_scheme = KratosSolid.ResidualBasedRotationNewmarkScheme(dynamic_factor, damp_factor_m)
            time_integration_method = KratosSolid.NewmarkStepMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            angular_time_integration_method = KratosSolid.NewmarkStepRotationMethod()
            angular_time_integration_method.AddToProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, angular_time_integration_method, self.process_info)
            angular_time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementRotationNewmarkScheme()
        elif(integration_method == "RotationBossak"):
            bossak_factor = self.implicit_solver_settings["bossak_factor"].GetDouble()
            self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = bossak_factor;
            #dynamic_factor = self.implicit_solver_settings["dynamic_factor"].GetDouble() # 0,1
            #mechanical_scheme = KratosSolid.ResidualBasedRotationNewmarkScheme(dynamic_factor, bossak_factor)
            time_integration_method = KratosSolid.BossakStepMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            angular_time_integration_method = KratosSolid.BossakStepRotationMethod()
            angular_time_integration_method.AddToProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, angular_time_integration_method, self.process_info)
            angular_time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementRotationBossakScheme()                  
        elif(integration_method == "RotationSimo"):
            bossak_factor = self.implicit_solver_settings["bossak_factor"].GetDouble()
            self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = bossak_factor;
            #dynamic_factor = self.implicit_solver_settings["dynamic_factor"].GetDouble() # 0,1
            #mechanical_scheme = KratosSolid.ResidualBasedRotationSimoScheme(dynamic_factor, bossak_factor)
            time_integration_method = KratosSolid.SimoStepMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            angular_time_integration_method = KratosSolid.SimoStepRotationMethod()
            angular_time_integration_method.AddToProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, angular_time_integration_method, self.process_info)
            angular_time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementRotationSimoScheme()
        elif(integration_method == "RotationEMC"):
            #dynamic_factor = self.implicit_solver_settings["dynamic_factor"].GetDouble() # 0,1       
            #mechanical_scheme = KratosSolid.ResidualBasedRotationEMCScheme(dynamic_factor)
            time_integration_method = KratosSolid.EmcStepMethod()
            time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
            time_integration_method.SetParameters(self.process_info)
            angular_time_integration_method = KratosSolid.EmcStepRotationMethod()
            angular_time_integration_method.AddToProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, angular_time_integration_method, self.process_info)
            angular_time_integration_method.SetParameters(self.process_info)
            mechanical_scheme = KratosSolid.ResidualBasedDisplacementRotationEmcScheme()
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
        # KratosMultiphysics.LineSearchStrategy (alternative -> to test)
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

