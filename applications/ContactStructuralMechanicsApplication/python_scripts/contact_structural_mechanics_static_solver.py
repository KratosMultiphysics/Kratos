#import kratos core and applications
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import the implicit solver (the explicit one is derived from it)
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_static_solver

# Import auxiliary methods
from KratosMultiphysics.ContactStructuralMechanicsApplication import auxiliary_methods_solvers

# Import the contact convergence criteria factory
from KratosMultiphysics.ContactStructuralMechanicsApplication.contact_convergence_criteria_factory import ContactConvergenceCriteriaFactory

def CreateSolver(model, custom_settings):
    return ContactStaticMechanicalSolver(model, custom_settings)

class ContactStaticMechanicalSolver(structural_mechanics_static_solver.StaticMechanicalSolver):
    """The structural mechanics contact static solver.

    This class creates the mechanical solvers for contact static analysis. It currently
    supports line search, linear, arc-length, form-finding and Newton-Raphson
    strategies.

    Public member variables:
    arc_length_settings -- settings for the arc length method.

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
        auxiliary_methods_solvers.AuxiliarySetSettings(self.settings, self.contact_settings)

        # Setting echo level
        self.echo_level =  self.settings["echo_level"].GetInt()

        # Initialize the processes list
        self.processes_list = None

        # Initialize the post process
        self.post_process = None

        KM.Logger.PrintInfo("::[Contact Mechanical Static Solver]:: ", "Construction of ContactMechanicalSolver finished")

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        auxiliary_methods_solvers.AuxiliaryValidateSettings(self)

    def AddVariables(self):

        super().AddVariables()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliary_methods_solvers.AuxiliaryAddVariables(self.main_model_part, mortar_type)

        KM.Logger.PrintInfo("::[Contact Mechanical Static Solver]:: ", "Variables ADDED")

    def AddDofs(self):

        super().AddDofs()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliary_methods_solvers.AuxiliaryAddDofs(self.main_model_part, mortar_type)

        KM.Logger.PrintInfo("::[Contact Mechanical Static Solver]:: ", "DOF's ADDED")

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
        auxiliary_methods_solvers.AuxiliarySolve(mechanical_solution_strategy)

    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        return is_converged

    def ExecuteFinalizeSolutionStep(self):
        super().ExecuteFinalizeSolutionStep()
        if self.contact_settings["ensure_contact"].GetBool():
            computing_model_part = self.GetComputingModelPart()
            CSMA.ContactUtilities.CheckActivity(computing_model_part)

    def ComputeDeltaTime(self):
        return auxiliary_methods_solvers.AuxiliaryComputeDeltaTime(self.main_model_part, self.GetComputingModelPart(), self.settings, self.contact_settings)

    def AddProcessesList(self, processes_list):
        self.processes_list = CSMA.ProcessFactoryUtility(processes_list)

    def AddPostProcess(self, post_process):
        self.post_process = CSMA.ProcessFactoryUtility(post_process)

    #### Private functions ####

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        return auxiliary_methods_solvers.AuxiliaryCreateConvergenceParameters(self.main_model_part, self.settings, self.contact_settings)

    def _CreateConvergenceCriterion(self):
        convergence_criterion = ContactConvergenceCriteriaFactory(self.main_model_part, self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _CreateLinearSolver(self):
        linear_solver = super()._CreateLinearSolver()
        return auxiliary_methods_solvers.AuxiliaryCreateLinearSolver(self.main_model_part, self.settings, self.contact_settings, self.linear_solver_settings, linear_solver)

    def _CreateBuilderAndSolver(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            linear_solver = self._GetLinearSolver()
            builder_and_solver_type: str = self.settings["builder_and_solver_settings"]["type"].GetString()

            if builder_and_solver_type == "block":
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
            elif self.settings["analysis_type"].GetString() == "non_linear":
                # Create strategy
                if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
                    mechanical_solution_strategy = self._create_contact_newton_raphson_strategy()
                elif self.settings["solving_strategy_settings"]["type"].GetString() == "line_search":
                    mechanical_solution_strategy = self._create_contact_line_search_strategy()
                elif self.settings["solving_strategy_settings"]["type"].GetString() == "arc_length":
                    mechanical_solution_strategy = self._create_arc_length_strategy()
            else:
                err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
                err_msg += "Available options are: \"linear\", \"non_linear\""
                raise Exception(err_msg)
        else:
            mechanical_solution_strategy = super()._CreateSolutionStrategy()

        return mechanical_solution_strategy

    def _create_contact_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self._GetScheme()
        self.linear_solver = self._GetLinearSolver()
        self.mechanical_convergence_criterion = self._GetConvergenceCriterion()
        self.builder_and_solver = self._GetBuilderAndSolver()
        return auxiliary_methods_solvers.AuxiliaryLineSearch(computing_model_part, self.mechanical_scheme, self.linear_solver, self.mechanical_convergence_criterion, self.builder_and_solver, self.settings, self.contact_settings, self.processes_list, self.post_process)

    def _create_contact_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self._GetScheme()
        self.mechanical_convergence_criterion = self._GetConvergenceCriterion()
        self.builder_and_solver = self._GetBuilderAndSolver()
        return auxiliary_methods_solvers.AuxiliaryNewton(computing_model_part, self.mechanical_scheme, self.mechanical_convergence_criterion, self.builder_and_solver, self.settings, self.contact_settings, self.processes_list, self.post_process)

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = auxiliary_methods_solvers.AuxiliaryContactSettings()
        this_defaults.RecursivelyAddMissingParameters(super(ContactStaticMechanicalSolver, cls).GetDefaultParameters())
        return this_defaults
