from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)


class ImplicitMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics implicit dynamic solver.

    This class creates the mechanical solvers for implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_type"   : "newmark",
            "damp_factor_m" :-0.3,
            "rayleigh_alpha": 0.0,
            "rayleigh_beta" : 0.0
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)
        print("::[ImplicitMechanicalSolver]:: Construction finished")

    def AddVariables(self):
        super(ImplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        print("::[ImplicitMechanicalSolver]:: Variables ADDED")
    
    def AddDofs(self):
        super(ImplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        print("::[ImplicitMechanicalSolver]:: DOF's ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.dynamic_settings["scheme_type"].GetString()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.dynamic_settings["rayleigh_beta"].GetDouble()
        if(scheme_type == "newmark"):
            damp_factor_m = 0.0
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif(scheme_type == "bossak"):
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif(scheme_type == "relaxation"):
            damp_factor_f =-0.3
            dynamic_factor_m = 10.0
            mechanical_scheme = StructuralMechanicsApplication.ResidualBasedRelaxationScheme(
                                                                       damp_factor_f, dynamic_factor_m)
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\", \"relaxation\""
            raise Exception(err_msg)
        return mechanical_scheme
