from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return ExplicitMechanicalSolver(main_model_part, custom_settings)

class ExplicitMechanicalSolver(BaseSolver.MechanicalSolver):
    """The solid mechanics explicit dynamic solver.

    This class creates the mechanical solvers for explicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See solid_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings): 

         # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "time_step_prediction_level": 0, 
            "max_delta_time": 1.0e-5, 
            "fraction_delta_time": 0.9, 
            "rayleigh_damping": false, 
            "rayleigh_alpha": 0.0,
            "rayleigh_beta" : 0.0

        }
        """)
        self._validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("CentralDifferences")
        
        # Construct the base solver.
        super(ExplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)

        print("::[Explicit Dynamics Solver]:: Constructed")       
   
    def SetVariables(self):

        BaseSolver.MechanicalSolver.SetVariables(self)

        integration_method = self.settings["time_integration_method"].GetString()
        # Add specific variables for the explicit time integration scheme
        if(integration_method == "Explicit"):
            self.nodal_variables = self.nodal_variables + ['NODAL_MASS','MIDDLE_VELOCITY','FORCE_RESIDUAL']
            # Add specific variables for the explicit time integration scheme in rotations
            if(self._check_input_dof("ROTATION") == True):
                self.nodal_variables = self.nodal_variables + ['INERTIA_DYADIC','MOMENT_RESIDUAL','POSITION_MOMENTUM','ROTATION_MOMENTUM', 'RESIDUAL_LYAPUNOV', 'TANGENT_LYAPUNOV']
                
        print("::[Explicit_Mechanical_Solver]:: Explicit Variables ADDED")


    #### Specific internal functions ####
    def _create_solution_scheme(self):

        scheme_type = self.settings["scheme_type"].GetString()
        
        if(scheme_type == "CentralDifferences"):
            mechanical_scheme = KratosSolid.ExplicitCentralDifferencesScheme(self.dynamic_settings["max_delta_time"].GetDouble(), 
                                                                             self.dynamic_settings["fraction_delta_time"].GetDouble(), 
                                                                             self.dynamic_settings["time_step_prediction_level"].GetDouble(), 
                                                                             self.dynamic_settings["rayleigh_damping"].GetBool())        
        else:
            raise Exception(" Scheme Type:", self.settings["scheme_type"].GetString()," not implemented yet.")
          
        return mechanical_scheme


    def _create_mechanical_solver(self):

        mechanical_solver = self._create_explicit_strategy()
        
        mechanical_solver.SetRebuildLevel(0) # 1 to recompute the mass matrix in each explicit step
        
        return mechanical_solver

    
    def _create_explicit_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        return KratosSolid.ExplicitStrategy(computing_model_part,
                                            mechanical_scheme, 
                                            linear_solver, 
                                            self.settings["compute_reactions"].GetBool(), 
                                            self.settings["reform_dofs_at_each_step"].GetBool(), 
                                            self.settings["move_mesh_flag"].GetBool())

 

