#import kratos core and applications
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import the implicit solver (the explicit one is derived from it)
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_implicit_dynamic_solver

# Import auxiliar methods
from KratosMultiphysics.ContactStructuralMechanicsApplication import auxiliar_methods_solvers

# Import the contact convergence criteria factory
from KratosMultiphysics.ContactStructuralMechanicsApplication.contact_convergence_criteria_factory import ContactConvergenceCriteriaFactory

def CreateSolver(model, custom_settings):
    return ContactImplicitMechanicalSolver(model, custom_settings)

class ContactImplicitMechanicalSolver(structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    """The structural mechanics contact implicit dynamic solver.

    This class creates the mechanical solvers for contact implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually

        # Construct the base solver.
        super().__init__(model, custom_settings)

        self.contact_settings = self.settings["contact_settings"]

        # Linear solver settings
        if self.settings.Has("linear_solver_settings"):
            self.linear_solver_settings = self.settings["linear_solver_settings"]
        else:
            self.linear_solver_settings = KM.Parameters("""{}""")

        # Setting default configurations true by default
        auxiliar_methods_solvers.AuxiliarSetSettings(self.settings, self.contact_settings)

        # Setting echo level
        self.echo_level =  self.settings["echo_level"].GetInt()

        # Initialize the processes list
        self.processes_list = None

        # Initialize the post process
        self.post_process = None

        KM.Logger.PrintInfo("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Construction of ContactMechanicalSolver finished")

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        auxiliar_methods_solvers.AuxiliarValidateSettings(self)

    def AddVariables(self):

        super().AddVariables()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliar_methods_solvers.AuxiliarAddVariables(self.main_model_part, mortar_type)

        KM.Logger.PrintInfo("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Variables ADDED")

    def AddDofs(self):

        super().AddDofs()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliar_methods_solvers.AuxiliarAddDofs(self.main_model_part, mortar_type)

        KM.Logger.PrintInfo("::[Contact Mechanical Implicit Dynamic Solver]:: ", "DOF's ADDED")

    def Initialize(self):
        super().Initialize() # The mechanical solver is created here.

        # No verbosity from strategy
        if self.contact_settings["silent_strategy"].GetBool():
            mechanical_solution_strategy = self._GetSolutionStrategy()
            mechanical_solution_strategy.SetEchoLevel(0)

        # We set the flag INTERACTION
        computing_model_part = self.GetComputingModelPart()
        if self.contact_settings["simplified_semi_smooth_newton"].GetBool():
            computing_model_part.Set(KM.INTERACTION, False)
        else:
            computing_model_part.Set(KM.INTERACTION, True)

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        mechanical_solution_strategy = self._GetSolutionStrategy()
        auxiliar_methods_solvers.AuxiliarSolve(mechanical_solution_strategy)

    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        return is_converged

    def ExecuteFinalizeSolutionStep(self):
        super().ExecuteFinalizeSolutionStep()
        if self.contact_settings["ensure_contact"].GetBool():
            computing_model_part = self.GetComputingModelPart()
            CSMA.ContactUtilities.CheckActivity(computing_model_part)

    def ComputeDeltaTime(self):
        return auxiliar_methods_solvers.AuxiliarComputeDeltaTime(self.main_model_part, self.GetComputingModelPart(), self.settings, self.contact_settings)

    def AddProcessesList(self, processes_list):
        self.processes_list = CSMA.ProcessFactoryUtility(processes_list)

    def AddPostProcess(self, post_process):
        self.post_process = CSMA.ProcessFactoryUtility(post_process)

    #### Private functions ####

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        return auxiliar_methods_solvers.AuxiliarCreateConvergenceParameters(self.main_model_part, self.settings, self.contact_settings)

    def _CreateConvergenceCriterion(self):
        convergence_criterion = ContactConvergenceCriteriaFactory(self.main_model_part, self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _CreateLinearSolver(self):
        linear_solver = super()._CreateLinearSolver()
        return auxiliar_methods_solvers.AuxiliarCreateLinearSolver(self.main_model_part, self.settings, self.contact_settings, self.linear_solver_settings, linear_solver)

    def _CreateBuilderAndSolver(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            linear_solver = self._GetLinearSolver()
            if self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
                builder_and_solver = CSMA.ContactResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                    # We use the elimination builder and solver
                    if self.settings["multi_point_constraints_used"].GetBool():
                        if (self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0):
                            self.GetComputingModelPart().Set(KM.TO_SPLIT) # We set the flag for some operations
                        builder_and_solver = CSMA.ContactResidualBasedEliminationBuilderAndSolverWithConstraints(linear_solver)
                    else:
                        builder_and_solver = CSMA.ContactResidualBasedEliminationBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = super()._CreateBuilderAndSolver()

        return builder_and_solver

    def _CreateSolutionStrategy(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            if self.settings["analysis_type"].GetString() == "linear":
                mechanical_solution_strategy = self._create_linear_strategy()
            else:
                if(self.settings["line_search"].GetBool()):
                    mechanical_solution_strategy = self._create_contact_line_search_strategy()
                else:
                    mechanical_solution_strategy = self._create_contact_newton_raphson_strategy()
        else:
            mechanical_solution_strategy = super()._CreateSolutionStrategy()

        return mechanical_solution_strategy

    def _create_contact_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self._GetScheme()
        self.linear_solver = self._GetLinearSolver()
        self.mechanical_convergence_criterion = self._GetConvergenceCriterion()
        self.builder_and_solver = self._GetBuilderAndSolver()
        return auxiliar_methods_solvers.AuxiliarLineSearch(computing_model_part, self.mechanical_scheme, self.linear_solver, self.mechanical_convergence_criterion, self.builder_and_solver, self.settings, self.contact_settings, self.processes_list, self.post_process)

    def _create_contact_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self._GetScheme()
        self.mechanical_convergence_criterion = self._GetConvergenceCriterion()
        self.builder_and_solver = self._GetBuilderAndSolver()
        return auxiliar_methods_solvers.AuxiliarNewton(computing_model_part, self.mechanical_scheme, self.mechanical_convergence_criterion, self.builder_and_solver, self.settings, self.contact_settings, self.processes_list, self.post_process)

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = auxiliar_methods_solvers.AuxiliarContactSettings()
        this_defaults.RecursivelyAddMissingParameters(super(ContactImplicitMechanicalSolver, cls).GetDefaultParameters())
        return this_defaults
