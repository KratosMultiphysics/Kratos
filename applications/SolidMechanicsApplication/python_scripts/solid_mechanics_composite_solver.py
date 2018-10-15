from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_segregated_solver as BaseSolver

def CreateSolver(custom_settings, Model):
    return CompositeSolver(Model, custom_settings)

#Base class to develop other solvers
class CompositeSolver(BaseSolver.SegregatedSolver):
    """The solid mechanics composite solver

    This class creates the a list of segregated solvers for the configuration of a composite strategy

    See solid_mechanics_segregated_solver.py for more information.
    """
    def __init__(self, Model, custom_settings):

        super(CompositeSolver, self).__init__(Model, custom_settings)

        # Composite solver counter
        self.solver_counter = 0

    def Clear(self):
        for solver in self.solvers:
            solver.Clear()

    def Check(self):
        for solver in self.solvers:
            solver.Check()

    #### Solve loop methods ####

    def Solve(self):

        # initialize solution step for all solvers
        for solver in self.solvers:
            solver.InitializeSolutionStep()

        # solver solution step for all solvers
        is_converged = True
        for solver in self.solvers:
            if not solver.SolveSolutionStep():
                is_converged = False
                break

        if is_converged is True:
            # finalize solution step for all solvers
            # allows to apply implex and other coupled updates
            for solver in self.solvers:
                solver.FinalizeSolutionStep()

        return is_converged

    # step by step:

    def InitializeSolutionStep(self):
        self.solvers[self.solver_counter].InitializeSolutionStep()

    def SolveSolutionStep(self):
        is_converged  = self.solvers[self.solver_counter].SolveSolutionStep()

        if(self.solver_counter == len(self.solvers)):
            self.solver_counter = 0
        else:
            self.solver_counter += 1

        return is_converged

    def FinalizeSolutionStep(self):
        self.solvers[self.solver_counter].FinalizeSolutionStep()

        if(self.solver_counter == len(self.solvers)):
            self.solver_counter = 0
        else:
            self.solver_counter += 1


    #### Solver internal methods ####

    #
    def _set_model_info(self):

        # Get solving model part from one of the solvers
        first_solver = self.settings["solvers"][0]["Parameters"]
        self.model_part = self.model[first_solver["solving_model_part"].GetString()]

        # Main model part from computing model part
        self.main_model_part = self.model_part.GetRootModelPart()

        # Process information
        self.process_info = self.main_model_part.ProcessInfo
