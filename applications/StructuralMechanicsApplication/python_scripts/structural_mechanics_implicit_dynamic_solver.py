from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(model, custom_settings):
    return ImplicitMechanicalSolver(model, custom_settings)


class ImplicitMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics implicit dynamic solver.

    This class creates the mechanical solvers for implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_type"   : "bossak",
            "damp_factor_m" :-0.3,
            "rayleigh_alpha": 0.0,
            "rayleigh_beta" : 0.0
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[ImplicitMechanicalSolver]:: ", "Construction finished")

        # Setting minimum buffer
        scheme_type = self.dynamic_settings["scheme_type"].GetString()
        if("bdf" in scheme_type or scheme_type == "backward_euler"):
            order = self._bdf_integration_order()
            self.settings["buffer_size"].SetInt(order + 1)

    def AddVariables(self):
        super(ImplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        self.print_on_rank_zero("::[ImplicitMechanicalSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        super(ImplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        self.print_on_rank_zero("::[ImplicitMechanicalSolver]:: ", "DOF's ADDED")

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
        elif(scheme_type == "pseudo_static"):
            mechanical_scheme = KratosMultiphysics.ResidualBasedPseudoStaticDisplacementScheme(StructuralMechanicsApplication.RAYLEIGH_BETA)
        elif(scheme_type.startswith("bdf") or scheme_type == "backward_euler"):
            order = self._bdf_integration_order()
            # In case of rotation dof we declare the dynamic variables
            if self.settings["rotation_dofs"].GetBool():
                dynamic_variables = KratosMultiphysics.Parameters(""" {
                    "variable"              : ["DISPLACEMENT","ROTATION"],
                    "first_derivative"      : ["VELOCITY","ANGULAR_VELOCITY"],
                    "second_derivative"     : ["ACCELERATION","ANGULAR_ACCELERATION"]
                    } """)
                mechanical_scheme = KratosMultiphysics.ResidualBasedBDFCustomScheme(order, dynamic_variables)
            else:
                mechanical_scheme = KratosMultiphysics.ResidualBasedBDFDisplacementScheme(order)
        elif(scheme_type == "relaxation"):
            damp_factor_f =-0.3
            dynamic_factor_m = 10.0
            mechanical_scheme = StructuralMechanicsApplication.ResidualBasedRelaxationScheme(
                                                                       damp_factor_f, dynamic_factor_m)
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\", \"pseudo_static\", \"backward_euler\", \"bdf1\", \"bdf2\", \"bdf3\", \"bdf4\", \"bdf5\", \"relaxation\""
            raise Exception(err_msg)
        return mechanical_scheme

    def _bdf_integration_order(self):
        scheme_type = self.dynamic_settings["scheme_type"].GetString()
        if (scheme_type == "backward_euler"):
            order = 1
        else:
            # BDF schemes can be from 1 to 5 order, so in order to detect the integration order from the scheme_type we remove the "bdf" string, that is, if the user tells bdf3 only 3 will remain when we remove bdf which corresponds to the method of choice
            order = int(scheme_type.replace("bdf", ""))

        # Warning
        if (order > 2):
            KratosMultiphysics.Logger.PrintWarning("WARNING:: BDF Order: ", str(order) + " constant time step must be considered")

        return order
