from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_segregated_solver as BaseSolver

def CreateSolver(custom_settings):
    return CompositeSolver(custom_settings)

#Base class to develop other solvers
class CompositeSolver(BaseSolver.SegregatedSolver):
    """The solid mechanics composite solver

    This class creates the a list of segregated solvers for the configuration of a composite strategy

    See solid_mechanics_segregated_solver.py for more information.
    """
    def __init__(self, custom_settings):

        super(CompositeSolver, self).__init__(custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solvers":[],
            "computing_parts": [],
            "processes":[]
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Create solvers list
        self.solvers = []
        solvers_list = self.settings["solvers"]
        for i in range(solvers_list.size()):
            solver_module = __import__(solvers_list[i]["solver_type"].GetString())
            self.solvers.append(solver_module.CreateSolver(solvers_list[i]["Parameters"]))

        # Computing parts (must be defined for each solver)
        if(solvers_list.size() != computing_parts.size() ):
            raise Exception( "Computing parts and solvers list must have the same size in a Composite Solver" )

        # Composite solver counter
        self.solver_counter = 0

        # Echo level
        self.echo_level = 0

        # Solver processes
        self.processes = []


    def SetComputingModelPart(self, computing_model_part):
        self.model_part = computing_model_part

        for create_part in self.create_parts:
            solver.SetComputingModelPart(self.model_part.GetSubModelPart(create_part["Parameters"]["model_part_name"].GetString()))

    def ExecuteInitialize(self):
        self._create_computing_sub_parts()
        super(CompositeSolver, self).ExecuteInitialize()

    #### Solve loop methods ####

    def Solve(self):
        self._create_computing_sub_parts()
        for solver in self.solvers:
            solver.Solve()

    # step by step:

    def InitializeSolutionStep(self):
        self._create_computing_sub_parts()
        self.solvers[self.solver_counter].InitializeSolutionStep()

    def SolveSolutionStep(self):
        self.solvers[self.solver_counter].SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.solvers[self.solver_counter].FinalizeSolutionStep()

        if(self.solver_counter == len(self.solvers)):
            self.solver_counter = 0
        else:
            self.solver_counter += 1


    #### Solver internal methods ####

    #
    def _create_computing_parts_process(self):
        for i in range(0,computing_parts.size()):
            create_computing_part = KratosSolid.TransferComputingModelPartProcess(self.model_part,computing_parts["Parameters"])
            self.create_parts.append(create_computing_part)

    #
    def _create_computing_sub_parts(self):
        for create_part in self.create_parts:
            create_part.Execute()
