from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
import trilinos_structural_mechanics_solver


def CreateSolver(model, custom_settings):
    return TrilinosImplicitMechanicalSolver(model, custom_settings)


class TrilinosImplicitMechanicalSolver(trilinos_structural_mechanics_solver.TrilinosMechanicalSolver):
    """The trilinos structural mechanics implicit dynamic solver.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    For more information see:
    structural_mechanics_solver.py
    trilinos_structural_mechanics_solver.py
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
        super(TrilinosImplicitMechanicalSolver, self).__init__(model, custom_settings)

    def AddVariables(self):
        super(TrilinosImplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        self.print_on_rank_zero("::[TrilinosImplicitMechanicalSolver]:: Variables ADDED")

    def AddDofs(self):
        super(TrilinosImplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        self.print_on_rank_zero("::[TrilinosImplicitMechanicalSolver]:: DOF's ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.dynamic_settings["scheme_type"].GetString()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.dynamic_settings["rayleigh_beta"].GetDouble()
        if (scheme_type == "newmark"):
            damp_factor_m = 0.0
        elif (scheme_type == "bossak"):
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\""
            raise Exception(err_msg)
        mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(damp_factor_m)
        return mechanical_scheme
