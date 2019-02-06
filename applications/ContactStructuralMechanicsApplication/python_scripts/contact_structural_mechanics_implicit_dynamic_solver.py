from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_implicit_dynamic_solver

# Import auxiliar methods
import auxiliar_methods_solvers

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

        ## Settings string in json format
        contact_settings = auxiliar_methods_solvers.AuxiliarContactSettings()

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.validate_and_transfer_matching_settings(self.settings, contact_settings)
        self.contact_settings = contact_settings["contact_settings"]

        # Linear solver settings
        if self.settings.Has("linear_solver_settings"):
            self.linear_solver_settings = self.settings["linear_solver_settings"]
        else:
            self.linear_solver_settings = KM.Parameters("""{}""")

        # Construct the base solver.
        super(ContactImplicitMechanicalSolver, self).__init__(model, self.settings)

        # Setting default configurations true by default
        auxiliar_methods_solvers.AuxiliarSetSettings(self.settings, self.contact_settings)

        # Setting echo level
        self.echo_level =  self.settings["echo_level"].GetInt()

        # Initialize the processes list
        self.processes_list = None

        # Initialize the post process
        self.post_process = None

        self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Construction of ContactMechanicalSolver finished")

    def AddVariables(self):

        super(ContactImplicitMechanicalSolver, self).AddVariables()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliar_methods_solvers.AuxiliarAddVariables(self.main_model_part, mortar_type)

        self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Variables ADDED")

    def AddDofs(self):

        super(ContactImplicitMechanicalSolver, self).AddDofs()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        auxiliar_methods_solvers.AuxiliarAddDofs(self.main_model_part, mortar_type)

        self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "DOF's ADDED")

    def Initialize(self):
        super(ContactImplicitMechanicalSolver, self).Initialize() # The mechanical solver is created here.

        # No verbosity from strategy
        if self.contact_settings["silent_strategy"].GetBool() is True:
            mechanical_solution_strategy = self.get_mechanical_solution_strategy()
            mechanical_solution_strategy.SetEchoLevel(0)

        # We set the flag INTERACTION
        computing_model_part = self.GetComputingModelPart()
        if self.contact_settings["simplified_semi_smooth_newton"].GetBool() is True:
            computing_model_part.Set(KM.INTERACTION, False)
        else:
            computing_model_part.Set(KM.INTERACTION, True)

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        mechanical_solution_strategy = self.get_mechanical_solution_strategy()
        auxiliar_methods_solvers.AuxiliarSolve(mechanical_solution_strategy)

    def SolveSolutionStep(self):
        is_converged = self.get_mechanical_solution_strategy().SolveSolutionStep()
        return is_converged

    def ExecuteFinalizeSolutionStep(self):
        super(ContactImplicitMechanicalSolver, self).ExecuteFinalizeSolutionStep()
        if self.contact_settings["ensure_contact"].GetBool():
            computing_model_part = self.GetComputingModelPart()
            CSMA.ContactUtilities.CheckActivity(computing_model_part)

    def ComputeDeltaTime(self):
        delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        if self.contact_settings["inner_loop_adaptive"].GetBool():
            process_info = self.GetComputingModelPart().ProcessInfo
            if process_info.Has(CSMA.INNER_LOOP_ITERATION):
                inner_iterations = process_info[CSMA.INNER_LOOP_ITERATION]
                if inner_iterations > 1:
                    delta_time = delta_time/float(inner_iterations)
                    self.print_on_rank_zero("::[Contact Mechanical Static Solver]:: ", "Advancing with a reduced delta time of ", delta_time)
        return delta_time

    def AddProcessesList(self, processes_list):
        self.processes_list = CSMA.ProcessFactoryUtility(processes_list)

    def AddPostProcess(self, post_process):
        self.post_process = CSMA.ProcessFactoryUtility(post_process)

    def print_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KM.Logger.PrintInfo(" ".join(map(str,args)))

    def print_warning_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KM.Logger.PrintWarning(" ".join(map(str,args)))

    #### Private functions ####

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        return auxiliar_methods_solvers.AuxiliarCreateConvergenceParameters(self.main_model_part, self.settings, self.contact_settings)

    def _create_convergence_criterion(self):
        import contact_convergence_criteria_factory
        convergence_criterion = contact_convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _create_linear_solver(self):
        linear_solver = super(ContactImplicitMechanicalSolver, self)._create_linear_solver()
        return auxiliar_methods_solvers.AuxiliarCreateLinearSolver(self.main_model_part, self.settings, self.contact_settings, self.linear_solver_settings, linear_solver)

    def _create_builder_and_solver(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            linear_solver = self.get_linear_solver()
            if self.settings["block_builder"].GetBool():
                if self.settings["multi_point_constraints_used"].GetBool():
                    builder_and_solver = CSMA.ContactResidualBasedBlockBuilderAndSolverWithConstraints(linear_solver)
                else:
                    builder_and_solver = CSMA.ContactResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                builder_and_solver = super(ContactImplicitMechanicalSolver, self)._create_builder_and_solver()
        else:
            builder_and_solver = super(ContactImplicitMechanicalSolver, self)._create_builder_and_solver()

        return builder_and_solver

    def _create_mechanical_solution_strategy(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            if self.settings["analysis_type"].GetString() == "linear":
                mechanical_solution_strategy = self._create_linear_strategy()
            else:
                if(self.settings["line_search"].GetBool()):
                    mechanical_solution_strategy = self._create_contact_line_search_strategy()
                else:
                    mechanical_solution_strategy = self._create_contact_newton_raphson_strategy()
        else:
            mechanical_solution_strategy = super(ContactImplicitMechanicalSolver, self)._create_mechanical_solution_strategy()

        return mechanical_solution_strategy

    def _create_contact_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self.get_solution_scheme()
        self.linear_solver = self.get_linear_solver()
        self.mechanical_convergence_criterion = self.get_convergence_criterion()
        self.builder_and_solver = self.get_builder_and_solver()
        return auxiliar_methods_solvers.AuxiliarLineSearch(computing_model_part, self.mechanical_scheme, self.linear_solver, self.mechanical_convergence_criterion, self.builder_and_solver, self.settings, self.contact_settings, self.processes_list, self.post_process)

    def _create_contact_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self.get_solution_scheme()
        self.linear_solver = self.get_linear_solver()
        self.mechanical_convergence_criterion = self.get_convergence_criterion()
        self.builder_and_solver = self.get_builder_and_solver()
        return auxiliar_methods_solvers.AuxiliarNewton(computing_model_part, self.mechanical_scheme, self.linear_solver, self.mechanical_convergence_criterion, self.builder_and_solver, self.settings, self.contact_settings, self.processes_list, self.post_process)
