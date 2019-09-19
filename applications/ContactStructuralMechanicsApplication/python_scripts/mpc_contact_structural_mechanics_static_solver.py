from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_static_solver

# Import auxiliar methods
from KratosMultiphysics.ContactStructuralMechanicsApplication import auxiliar_methods_solvers

# Import convergence_criteria_factory
from KratosMultiphysics.StructuralMechanicsApplication import convergence_criteria_factory

def GetDefaults():
    this_defaults = KratosMultiphysics.Parameters("""
    {
        "mpc_contact_settings" :
        {
            "contact_type"                  : "Frictionless",
            "simplified_semi_smooth_newton" : false,
            "inner_loop_iterations"         : 10,
            "update_each_nl_iteration"      : false,
            "enforce_ntn"                   : false
        }
    }
    """)
    return this_defaults

def CreateSolver(model, custom_settings):
    return MPCContactStaticSolver(model, custom_settings)

class MPCContactStaticSolver(structural_mechanics_static_solver.StaticMechanicalSolver):
    """The MPC contact static solver.

    This class creates the mechanical solvers for contact static analysis. It currently
    supports linear and Newton-Raphson strategies.

    Public member variables:

    See structural_mechanics_solver.py for more information.
    """

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually

        # Construct the base solver.
        super(MPCContactStaticSolver, self).__init__(model, custom_settings)

        self.mpc_contact_settings = self.settings["mpc_contact_settings"]
        self.mpc_contact_settings.RecursivelyAddMissingParameters(GetDefaults()["mpc_contact_settings"])

        # Setting the parameters
        if not self.settings["compute_reactions"].GetBool():
            KratosMultiphysics.Logger.PrintInfo("Compute reactions", "Storage must be cleared each step. Switching to True")
            self.settings["compute_reactions"].SetBool(True)
        if not self.settings["clear_storage"].GetBool():
            KratosMultiphysics.Logger.PrintInfo("Clear storage", "Storage must be cleared each step. Switching to True")
            self.settings["clear_storage"].SetBool(True)
        if not self.settings["reform_dofs_at_each_step"].GetBool():
            KratosMultiphysics.Logger.PrintInfo("Reform DoFs", "DoF must be reformed each time step. Switching to True")
            self.settings["reform_dofs_at_each_step"].SetBool(True)

        # Construct the base solver.
        KratosMultiphysics.Logger.PrintInfo("::[MPCContactStaticSolver]:: ", "Construction finished")

    def AddVariables(self):

        super(MPCContactStaticSolver, self).AddVariables()

        # We add the contact related variables
        contact_type = self.mpc_contact_settings["contact_type"].GetString()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)  # Add normal
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H) # Add nodal size variable
        self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)  # Add normal contact gap
        if contact_type == "Frictional":
            self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_SLIP) # Add contact slip

    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("::[MPCContactStaticSolver]:: ", "Initializing ...")

        super(MPCContactStaticSolver, self).Initialize() # The mechanical solver is created here.

        # We set the flag INTERACTION
        if self.mpc_contact_settings["simplified_semi_smooth_newton"].GetBool():
            computing_model_part = self.GetComputingModelPart()
            computing_model_part.ProcessInfo.Set(KratosMultiphysics.INTERACTION, True)

        KratosMultiphysics.Logger.PrintInfo("::[MPCContactStaticSolver]:: ", "Finished initialization.")

    def ComputeDeltaTime(self):
        return auxiliar_methods_solvers.AuxiliarComputeDeltaTime(self.main_model_part, self.GetComputingModelPart(), self.settings, self.mpc_contact_settings)

    #### Private functions ####

    def _create_convergence_criterion(self):
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        conv_criteria = convergence_criterion.mechanical_convergence_criterion
        contact_criteria = ContactStructuralMechanicsApplication.MPCContactCriteria()
        return KratosMultiphysics.AndCriteria(conv_criteria, contact_criteria)

    def _create_mechanical_solution_strategy(self):
        mechanical_solution_strategy = self._create_contact_newton_raphson_strategy()
        return mechanical_solution_strategy

    def _create_contact_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self.get_solution_scheme()
        self.linear_solver = self.get_linear_solver()
        self.mechanical_convergence_criterion = self.get_convergence_criterion()
        self.builder_and_solver = self.get_builder_and_solver()
        newton_parameters = KratosMultiphysics.Parameters("""{}""")
        newton_parameters.AddValue("inner_loop_iterations", self.mpc_contact_settings["inner_loop_iterations"])
        newton_parameters.AddValue("update_each_nl_iteration", self.mpc_contact_settings["update_each_nl_iteration"])
        newton_parameters.AddValue("enforce_ntn", self.mpc_contact_settings["enforce_ntn"])
        return ContactStructuralMechanicsApplication.ResidualBasedNewtonRaphsonMPCContactStrategy(computing_model_part,
                                                                    self.mechanical_scheme,
                                                                    self.linear_solver,
                                                                    self.mechanical_convergence_criterion,
                                                                    self.builder_and_solver,
                                                                    self.settings["max_iteration"].GetInt(),
                                                                    self.settings["compute_reactions"].GetBool(),
                                                                    self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                    self.settings["move_mesh_flag"].GetBool(),
                                                                    newton_parameters
                                                                    )

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = GetDefaults()
        this_defaults.RecursivelyAddMissingParameters(super(MPCContactStaticSolver, cls).GetDefaultSettings())
        return this_defaults
