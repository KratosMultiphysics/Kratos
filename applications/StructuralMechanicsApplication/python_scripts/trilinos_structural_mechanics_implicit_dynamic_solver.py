from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication
import trilinos_structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return TrilinosImplicitMechanicalSolver(main_model_part, custom_settings)


class TrilinosImplicitMechanicalSolver(trilinos_structural_mechanics_solver.TrilinosMechanicalSolver):
    """The trilinos structural mechanics implicit dynamic solver.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    For more information see:
    structural_mechanics_solver.py
    trilinos_structural_mechanics_solver.py
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "damp_factor_m" :-0.3,
            "rayleigh_alpha": 0.0,
            "rayleigh_beta" : 0.0
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("newmark")
        # Construct the base solver.
        super(TrilinosImplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)

    def AddVariables(self):
        super(TrilinosImplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        print("::[TrilinosImplicitMechanicalSolver]:: Variables ADDED")
    
    def AddDofs(self):
        super(TrilinosImplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        print("::[TrilinosImplicitMechanicalSolver]:: DOF's ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_type"].GetString()
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
