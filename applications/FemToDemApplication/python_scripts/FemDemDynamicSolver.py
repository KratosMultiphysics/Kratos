from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.FemToDemApplication as KratosFemDem

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import FemDemMechanicalSolver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver(BaseSolver.FemDemMechanicalSolver):

    def __init__(self, main_model_part, custom_settings): 
        
        # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "damp_factor_m" :-0.01,
            "dynamic_factor": 1.0,
            "rayleigh_damping": false, 
            "rayleigh_alpha": 0.0,
            "rayleigh_beta" : 0.0
        }
        """)

        self._validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("Newmark")

        if not custom_settings.Has("extrapolation_required"):
            custom_settings.AddEmptyValue("extrapolation_required")
            custom_settings["extrapolation_required"].SetBool(False)

        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)

        print("::[Implicit_Dynamic_Solver]:: Constructed")

    #### Solver internal methods ####
    
    def _create_solution_scheme(self):
        
        scheme_type = self.settings["scheme_type"].GetString()
        
        if( self.dynamic_settings["rayleigh_damping"].GetBool() == True ):
            self.main_model_part.ProcessInfo[KratosSolid.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
            self.main_model_part.ProcessInfo[KratosSolid.RAYLEIGH_BETA]  = self.dynamic_settings["rayleigh_beta"].GetDouble()
        else:
            self.main_model_part.ProcessInfo[KratosSolid.RAYLEIGH_ALPHA] = 0.0
            self.main_model_part.ProcessInfo[KratosSolid.RAYLEIGH_BETA]  = 0.0
            
        if(self.settings["component_wise"].GetBool() == True):
            dynamic_factor = self.dynamic_settings["dynamic_factor"].GetDouble()        
            damp_factor_m  = self.dynamic_settings["damp_factor_m"].GetDouble()
            mechanical_scheme = KratosSolid.ComponentWiseBossakScheme(damp_factor_m, dynamic_factor)
        elif(scheme_type == "Newmark"):
            damp_factor_m = 0.0
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif(scheme_type == "Bossak"):
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif(scheme_type == "RotationNewmark"):
            dynamic_factor = self.dynamic_settings["dynamic_factor"].GetDouble() # 0,1 
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
            mechanical_scheme = KratosSolid.ResidualBasedRotationNewmarkScheme(dynamic_factor, damp_factor_m)
        elif(scheme_type == "RotationSimo"):
            dynamic_factor = self.dynamic_settings["dynamic_factor"].GetDouble() # 0,1       
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
            mechanical_scheme = KratosSolid.ResidualBasedRotationSimoScheme(dynamic_factor, damp_factor_m)
        elif(scheme_type == "RotationEMC"):
            dynamic_factor = self.dynamic_settings["dynamic_factor"].GetDouble() # 0,1       
            mechanical_scheme = KratosSolid.ResidualBasedRotationEMCScheme(dynamic_factor)
        else:
            raise Exception("Unsupported scheme_type: " + scheme_type)
                    
        return mechanical_scheme
    
    def _create_mechanical_solver(self):
        if(self.settings["component_wise"].GetBool() == True):
            mechanical_solver = self._create_component_wise_strategy()
        elif(self.settings["line_search"].GetBool() == True):
            if(self.settings["implex"].GetBool() == True):
                mechanical_solver = self._create_line_search_implex_strategy()
            else:
                mechanical_solver = self._create_line_search_strategy()
        else:
            if self.settings["extrapolation_required"].GetBool():
                mechanical_solver = self._create_newton_raphson_hexaedrons_strategy()
            else:
                mechanical_solver = self._create_newton_raphson_strategy()
        return mechanical_solver


    def _create_component_wise_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ComponentWiseNewtonRaphsonStrategy(computing_model_part, 
                                                              mechanical_scheme, 
                                                              linear_solver, 
                                                              mechanical_convergence_criterion, 
                                                              builder_and_solver, 
                                                              self.settings["max_iteration"].GetInt(),
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              self.settings["move_mesh_flag"].GetBool())
     
    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        # KratosMultiphysics.LineSearchStrategy (alternative -> to test)
        return KratosSolid.ResidualBasedNewtonRaphsonLineSearchStrategy(computing_model_part, 
                                                                        mechanical_scheme, 
                                                                        linear_solver, 
                                                                        mechanical_convergence_criterion, 
                                                                        builder_and_solver, 
                                                                        self.settings["max_iteration"].GetInt(),
                                                                        self.settings["compute_reactions"].GetBool(),
                                                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                        self.settings["move_mesh_flag"].GetBool())
    
    def _create_line_search_implex_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ResidualBasedNewtonRaphsonLineSearchImplexStrategy(computing_model_part, 
                                                                              mechanical_scheme, 
                                                                              linear_solver, 
                                                                              mechanical_convergence_criterion, 
                                                                              builder_and_solver, 
                                                                              self.settings["max_iteration"].GetInt(),
                                                                              self.settings["compute_reactions"].GetBool(),
                                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                              self.settings["move_mesh_flag"].GetBool())
    
    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part, 
                                                                     mechanical_scheme, 
                                                                     linear_solver, 
                                                                     mechanical_convergence_criterion, 
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_hexaedrons_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        return KratosFemDem.HexahedraNewtonRaphsonStrategy(computing_model_part, 
                                                                     mechanical_scheme, 
                                                                     linear_solver, 
                                                                     mechanical_convergence_criterion, 
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool()) 
                                        

