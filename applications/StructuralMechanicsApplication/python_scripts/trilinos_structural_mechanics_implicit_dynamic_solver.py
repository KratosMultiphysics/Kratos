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
            "damp_factor_m" :-0.01
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("Newmark")
        # Construct the base solver.
        super(TrilinosImplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)

    def AddVariables(self):
        super(TrilinosImplicitMechanicalSolver, self).AddVariables()
        # Add dynamic variables.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        print("::[TrilinosImplicitMechanicalSolver]:: Variables ADDED")

    #### Private functions ####

    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_type"].GetString()
        if (scheme_type == "Newmark"):
            damp_factor_m = 0.0
        elif (scheme_type == "Bossak"):
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
        else:
            raise Exception("Unsupported scheme_type: " + scheme_type)
        mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(damp_factor_m)
        return mechanical_scheme
