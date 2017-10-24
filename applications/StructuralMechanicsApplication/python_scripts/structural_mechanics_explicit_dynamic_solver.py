from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return ExplicitMechanicalSolver(main_model_part, custom_settings)

class ExplicitMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The solid mechanics explicit dynamic solver.

    This class creates the mechanical solvers for explicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the explicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
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
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("centraldifferences")
        
        # Construct the base solver.
        super(ExplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)
        print("::[ExplicitMechanicalSolver]:: Construction finished")    

    def AddVariables(self):
        super(ExplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        print("::[ExplicitMechanicalSolver]:: Variables ADDED")
    
    def AddDofs(self):
        super(ExplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.MIDDLE_VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.MIDDLE_VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.MIDDLE_VELOCITY_Z,self.main_model_part)
        print("::[ExplicitMechanicalSolver]:: DOF's ADDED")


    #### Specific internal functions ####
    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_type"].GetString()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.dynamic_settings["rayleigh_beta"].GetDouble()


        if(scheme_type == "centraldifferences"):
            mechanical_scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(self.dynamic_settings["max_delta_time"].GetDouble(), 
                                                                             self.dynamic_settings["fraction_delta_time"].GetDouble(), 
                                                                             self.dynamic_settings["time_step_prediction_level"].GetDouble(), 
                                                                             self.dynamic_settings["rayleigh_damping"].GetBool())                                                                     
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"centraldifferences\""
            raise Exception(err_msg)
        return mechanical_scheme


   
