import numpy as np

# Importing the Kratos Library
import KratosMultiphysics
# Import applications
import KratosMultiphysics.SPHApplication as SPH
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
# Import base class file
from KratosMultiphysics.SPHApplication.sph_solver import SPHSolver

def CreateSolver(model, custom_settings):
    return ImplicitSPHSolver(model, custom_settings)

class ImplicitSPHSolver(SPHSolver):
    
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitSPHSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "time_integration_method" : "implicit",
            "scheme_type"             : "newmark",
            "damp_factor_m"           : 0.0,
            "newmark_beta"            : 0.25,
            "rayleigh_alpha"          : 0.0,
            "rayleigh_beta"           : 0.0,
            "first_order_time_integrator_settings" : {
                "integration_variable_pairs" : []
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults
    
    def AddVariables(self):
        super().AddVariables()
        self._add_dynamic_variables()
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitSPHSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        super().AddDofs()
        self._add_dynamic_dofs()
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitSPHSolver]:: ", "DOF's ADDED")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        # Some pre-processes may affect the system of equations, we rebuild the equation ids
        process_info = self.main_model_part.ProcessInfo
        if process_info[KratosMultiphysics.STEP] == 1 and process_info[StructuralMechanicsApplication.RESET_EQUATION_IDS]:
            # Resetting the global equations ids
            self._GetBuilderAndSolver().SetUpSystem(self.GetComputingModelPart())
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.ExposeSystemMatrix()

    def _CreateScheme(self):
        scheme_type = self.settings["scheme_type"].GetString().strip().lower()

        # Setting the Rayleigh damping parameters
        process_info = self.main_model_part.ProcessInfo
        process_info[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.settings["rayleigh_alpha"].GetDouble()
        process_info[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.settings["rayleigh_beta"].GetDouble()

        # Setting the time integration schemes
        if scheme_type in ("newmark", "bossak"):
            scheme_settings = KratosMultiphysics.Parameters("""{
                "damp_factor_m" : 0.0,
                "newmark_beta" : 0.0,
                "projection_variables_list" : []
            }""")
            scheme_settings["damp_factor_m"].SetDouble(0.0 if scheme_type == "newmark" else self.settings["damp_factor_m"].GetDouble())
            scheme_settings["newmark_beta"].SetDouble(self.settings["newmark_beta"].GetDouble())
            sph_scheme = StructuralMechanicsApplication.StructuralMechanicsBossakScheme(scheme_settings)

#        elif scheme_type == "backward_euler":
#            if not self.settings["strain_dofs"].GetBool():
#                raise Exception('The "backward_euler" scheme is intended for the first-order mixed formulations and requires "strain_dofs" to be true.')
#
#            scheme_settings = self.settings["first_order_time_integrator_settings"].Clone()
#            scheme_settings.AddInt("domain_size", self.settings["domain_size"].GetInt())
#
#            integration_variable_pairs = scheme_settings["integration_variable_pairs"]
#            if integration_variable_pairs.size() == 0:
#                integration_variable_pairs.Append(KratosMultiphysics.Parameters("""{
#                    "primary_variable" : "VELOCITY",
#                    "first_derivative" : "ACCELERATION"
#                }"""))
#
#                deformation_gradient_components = ["XX", "YY", "XY", "YX"]
#                if self.settings["domain_size"].GetInt() == 3:
#                    deformation_gradient_components = ["XX", "YY", "ZZ", "XY", "XZ", "YX", "YZ", "ZX", "ZY"]
#
#                for component in deformation_gradient_components:
#                    pair_settings = KratosMultiphysics.Parameters("""{
#                        "primary_variable" : "",
#                        "first_derivative" : ""
#                    }""")
#                    pair_settings["primary_variable"].SetString(
#                        "DEFORMATION_GRADIENT_" + component)
#                    pair_settings["first_derivative"].SetString(
#                        "DEFORMATION_GRADIENT_DOT_" + component)
#                    integration_variable_pairs.Append(pair_settings)
#
#            sph_scheme = SPH.GenericBackwardEulerScheme(scheme_settings)
#
#        elif scheme_type == "crank_nicolson" : 
#            if not self.settings["strain_dofs"].GetBool():
#                raise Exception('The "crank_nicolson" scheme is intended for the first-order mixed formulations and requires "strain_dofs" to be true.')
#
#            scheme_settings = self.settings["first_order_time_integrator_settings"].Clone()
#            scheme_settings.AddInt("domain_size", self.settings["domain_size"].GetInt())
#
#            integration_variable_pairs = scheme_settings["integration_variable_pairs"]
#            if integration_variable_pairs.size() == 0:
#                integration_variable_pairs.Append(KratosMultiphysics.Parameters("""{
#                    "primary_variable" : "VELOCITY",
#                    "first_derivative" : "ACCELERATION"
#                }"""))
#
#                deformation_gradient_components = ["XX", "YY", "XY", "YX"]
#                if self.settings["domain_size"].GetInt() == 3:
#                    deformation_gradient_components = ["XX", "YY", "ZZ", "XY", "XZ", "YX", "YZ", "ZX", "ZY"]
#
#                for component in deformation_gradient_components:
#                    pair_settings = KratosMultiphysics.Parameters("""{
#                        "primary_variable" : "",
#                        "first_derivative" : ""
#                    }""")
#                    pair_settings["primary_variable"].SetString(
#                        "DEFORMATION_GRADIENT_" + component)
#                    pair_settings["first_derivative"].SetString(
#                        "DEFORMATION_GRADIENT_DOT_" + component)
#                    integration_variable_pairs.Append(pair_settings)
#
#            sph_scheme = SPH.GenericCrankNicolsonScheme(scheme_settings)
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\", \"pseudo_static\", \"backward_euler\""
            raise Exception(err_msg)
        return sph_scheme
