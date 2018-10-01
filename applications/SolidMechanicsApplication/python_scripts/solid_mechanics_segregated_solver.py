from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_monolithic_solver as BaseSolver

def CreateSolver(custom_settings, Model):
    return SegregatedSolver(Model, custom_settings)

#Base class to develop other solvers
class SegregatedSolver(BaseSolver.MonolithicSolver):
    """The solid mechanics segregated solver

    This class creates the a list of monolithic solvers for the configuration of a segregated strategy

    See solid_mechanics_monolithic_solver.py for more information.
    """
    def __init__(self, Model, custom_settings):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solving_model_part": "computing_domain",
            "solvers":[],
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
            self.solvers.append(solver_module.CreateSolver(solvers_list[i]["Parameters"], Model))

        # Model
        self.model = Model

        # Echo level
        self.echo_level = 0

        # Solver processes
        self.processes = []


    def ExecuteInitialize(self):
        for solver in self.solvers:
            solver._set_model_info()

        super(SegregatedSolver, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for solver in self.solvers:
            solver.ExecuteBeforeSolutionLoop()

        self.SetEchoLevel(self.echo_level)


    def GetVariables(self):
        nodal_variables = []
        for solver in self.solvers:
            nodal_variables = nodal_variables + solver.GetVariables()
        return nodal_variables

    def SetEchoLevel(self, level):
        for solver in self.solvers:
            solver.SetEchoLevel(level)
        self.echo_level = level

    def Clear(self):
        for solver in self.solvers:
            solver.Clear()

    #### Solver internal methods ####

    def _check_initialized(self):
        if( not self._is_not_restarted() ):
            for solver in self.solvers:
                if hasattr(mechanical_solver, 'SetInitializePerformedFlag'):
                    solver._get_mechanical_solver().SetInitializePerformedFlag(True)
                else:
                    solver._get_mechanical_solver().Set(KratosSolid.SolverLocalFlags.INITIALIZED, True)

            if hasattr(mechanical_solver, 'SetInitializePerformedFlag'):
                self._get_mechanical_solver().SetInitializePerformedFlag(True)
            else:
                self._get_mechanical_solver().Set(KratosSolid.SolverLocalFlags.INITIALIZED, True)

    def _create_mechanical_solver(self):
        strategies = []
        for solver in self.solvers:
            strategies.append(solver._get_mechanical_solver())
        options = KratosMultiphysics.Flags()
        mechanical_solver =  KratosSolid.SegregatedStrategy(self.model_part, options, strategies)
        mechanical_solver.Set(KratosSolid.SolverLocalFlags.ADAPTIVE_SOLUTION,self.settings["solving_strategy_settings"]["adaptive_solution"].GetBool())
        return mechanical_solver
    #
    def _get_dofs(self):
        dof_variables = []
        dof_reactions = []
        for solver in self.solvers:
            solver_dof_variables, solver_dof_reactions = solver._get_dofs()
            dof_variables = dof_variables + solver_dof_variables
            dof_reactions = dof_reactions + solver_dof_reactions

        return dof_variables, dof_reactions


    def _add_dofs(self):
        dof_variables, dof_reactions = self._get_dofs()
        AddDofsProcess = KratosSolid.AddDofsProcess(self.main_model_part, dof_variables, dof_reactions)
        AddDofsProcess.Execute()
        if( self.echo_level > 1 ):
            print(dof_variables + dof_reactions)
            print("::[-------Solver------]:: DOF's ADDED")

    #
    def _get_time_integration_methods(self):
        scalar_integration_methods = {}
        component_integration_methods = {}
        for solver in self.solvers:
            solver_scalar_integration_methods, solver_component_integration_methods = solver._get_time_integration_methods()
            scalar_integration_methods.update(solver_scalar_integration_methods)
            component_integration_methods.update(solver_component_integration_methods)

        return scalar_integration_methods, component_integration_methods

    #
    def _get_minimum_buffer_size(self):
        buffer_size = 2
        for solver in self.solvers:
            size = solver._get_minimum_buffer_size()
            if( size > buffer_size ):
                buffer_size = size
        return buffer_size;
