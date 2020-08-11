#import kratos core and applications
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import the explicit solver (the explicit one is derived from it)
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_explicit_dynamic_solver

# Import auxiliar methods
from KratosMultiphysics.ContactStructuralMechanicsApplication import auxiliar_methods_solvers

def CreateSolver(model, custom_settings):
    return ContactExplicitMechanicalSolver(model, custom_settings)

class ContactExplicitMechanicalSolver(structural_mechanics_explicit_dynamic_solver.ExplicitMechanicalSolver):
    """The structural mechanics contact explicit dynamic solver.

    This class creates the mechanical solvers for contact explicit dynamic analysis.
    It currently supports central difference method

    Public member variables:
    dynamic_settings -- settings for the explicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(ContactExplicitMechanicalSolver, self).__init__(model, custom_settings)

        self.contact_settings = self.settings["contact_settings"]

        # Setting default configurations true by default
        auxiliar_methods_solvers.AuxiliarSetSettings(self.settings, self.contact_settings)

        # Getting delta_time_factor_for_contact
        self.delta_time_factor_for_contact = self.contact_settings["delta_time_factor_for_contact"].GetDouble()

        # Setting echo level
        self.echo_level =  self.settings["echo_level"].GetInt()

        KM.Logger.PrintInfo("::[Contact Mechanical Explicit Dynamic Solver]:: ", "Construction of ContactMechanicalSolver finished")

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        auxiliar_methods_solvers.AuxiliarValidateSettings(self)

    def AddVariables(self):

        super(ContactExplicitMechanicalSolver, self).AddVariables()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliar_methods_solvers.AuxiliarAddVariables(self.main_model_part, mortar_type)

        KM.Logger.PrintInfo("::[Contact Mechanical Explicit Dynamic Solver]:: ", "Variables ADDED")

    def AddDofs(self):

        super(ContactExplicitMechanicalSolver, self).AddDofs()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliar_methods_solvers.AuxiliarAddDofs(self.main_model_part, mortar_type)

        KM.Logger.PrintInfo("::[Contact Mechanical Explicit Dynamic Solver]:: ", "DOF's ADDED")

    def Initialize(self):
        super(ContactExplicitMechanicalSolver, self).Initialize() # The mechanical solver is created here.

        # No verbosity from strategy
        if self.contact_settings["silent_strategy"].GetBool():
            mechanical_solution_strategy = self.get_mechanical_solution_strategy()
            mechanical_solution_strategy.SetEchoLevel(0)

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        mechanical_solution_strategy = self.get_mechanical_solution_strategy()
        auxiliar_methods_solvers.AuxiliarSolve(mechanical_solution_strategy)

    def SolveSolutionStep(self):
        is_converged = self.get_mechanical_solution_strategy().SolveSolutionStep()
        return is_converged

    def ExecuteFinalizeSolutionStep(self):
        super(ContactExplicitMechanicalSolver, self).ExecuteFinalizeSolutionStep()
        if self.contact_settings["ensure_contact"].GetBool():
            computing_model_part = self.GetComputingModelPart()
            CSMA.ContactUtilities.CheckActivity(computing_model_part)

    def ComputeDeltaTime(self):
        self.delta_time = super(ContactExplicitMechanicalSolver, self).ComputeDeltaTime()
        if self.GetComputingModelPart().Is(KM.CONTACT):
            return self.delta_time_factor_for_contact * self.delta_time
        else:
            return self.delta_time

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = auxiliar_methods_solvers.AuxiliarExplicitContactSettings()
        this_defaults.RecursivelyAddMissingParameters(super(ContactExplicitMechanicalSolver, cls).GetDefaultSettings())
        return this_defaults
