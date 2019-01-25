from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver

# Import common utilities
from . import common_methods_solvers

def CreateSolver(model, custom_settings):
    return ExplicitMechanicalSolver(model, custom_settings)

class ExplicitMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics explicit dynamic solver.

    This class creates the mechanical solvers for explicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the explicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
         # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_type"                : "central_differences",
            "time_step_prediction_level" : 0,
            "max_delta_time"             : 1.0e-5,
            "fraction_delta_time"        : 0.9,
            "rayleigh_alpha"             : 0.0,
            "rayleigh_beta"              : 0.0
        }
        """)
        damping_settings = common_methods_solvers.AuxiliarDampingSettings()
        if custom_settings.Has("damping_settings"):
            if custom_settings["damping_settings"].Has("determine_rayleigh_damping_settings"):
                if custom_settings["damping_settings"]["determine_rayleigh_damping_settings"].Has("eigen_system_settings"):
                    if custom_settings["damping_settings"]["determine_rayleigh_damping_settings"]["eigen_system_settings"].Has("solver_type"):
                        solver_type = custom_settings["damping_settings"]["determine_rayleigh_damping_settings"]["eigen_system_settings"]["solver_type"].GetString()
                        eigen_system_settings = common_methods_solvers.AuxiliarEigenSettings(solver_type)
                        damping_settings["damping_settings"]["determine_rayleigh_damping_settings"]["eigen_system_settings"] = eigen_system_settings["eigen_system_settings"]
        self.dynamic_settings.AddValue("damping_settings", damping_settings["damping_settings"])
        self.validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(ExplicitMechanicalSolver, self).__init__(model, custom_settings)
        # Lumped mass-matrix is necessary for explicit analysis
        self.main_model_part.ProcessInfo[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = True
        self.print_on_rank_zero("::[ExplicitMechanicalSolver]:: Construction finished")

    def AddVariables(self):
        super(ExplicitMechanicalSolver, self).AddVariables()
        self._add_dynamic_variables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)

        if (self.settings["rotation_dofs"].GetBool()):
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)

        self.print_on_rank_zero("::[ExplicitMechanicalSolver]:: Variables ADDED")

    def AddDofs(self):
        super(ExplicitMechanicalSolver, self).AddDofs()
        self._add_dynamic_dofs()
        self.print_on_rank_zero("::[ExplicitMechanicalSolver]:: DOF's ADDED")

    def InitializeSolutionStep(self):
        # Using the base InitializeSolutionStep
        super(ExplicitMechanicalSolver, self).InitializeSolutionStep()

        # We compute for the different parts (each body can have it own damping)
        if self.dynamic_settings["damping_settings"]["determine_rayleigh_damping"].GetBool() and self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] == 1:
            common_methods_solvers.ComputeDampingCoefficients(self.model, self.settings, self.dynamic_settings["damping_settings"])

    #### Specific internal functions ####
    def _create_solution_scheme(self):
        scheme_type = self.dynamic_settings["scheme_type"].GetString()

        # Setting the Rayleigh damping parameters
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
        self.main_model_part.ProcessInfo[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.dynamic_settings["rayleigh_beta"].GetDouble()

        # Setting the time integration schemes
        if(scheme_type == "central_differences"):
            mechanical_scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(self.dynamic_settings["max_delta_time"].GetDouble(),
                                                                             self.dynamic_settings["fraction_delta_time"].GetDouble(),
                                                                             self.dynamic_settings["time_step_prediction_level"].GetDouble())
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"central_differences\""
            raise Exception(err_msg)
        return mechanical_scheme

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()

        mechanical_solution_strategy = StructuralMechanicsApplication.MechanicalExplicitStrategy(computing_model_part,
                                            mechanical_scheme,
                                            self.settings["compute_reactions"].GetBool(),
                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                            self.settings["move_mesh_flag"].GetBool())

        mechanical_solution_strategy.SetRebuildLevel(0)
        return mechanical_solution_strategy

    #### Private functions ####

