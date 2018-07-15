from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
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

        # Composite solver counter
        self.solver_counter = 0

        # Solver processes
        self.processes = []
        for i in range(0,process_list.size()):
            self.processes.append(process_list[i])

        
    def SetComputingModelPart(self, computing_model_part):
        self.model_part = computing_model_part

        counter = 0
        for solver in self.solvers:
            solver.SetComputingModelPart(solver._create_computing_sub_model_part(computing_model_part,counter))
            counter+=1

    def ExecuteInitialize(self):
        super(CompositeSolver, self).ExecuteInitialize()
        self._processes_execute_initialize()
        
    #### Solve loop methods ####

    def Solve(self):
        for solver in self.solvers:
            self._processes_execute_initialize_solution_step()
            solver.Solve()
            self._processes_execute_finalize_solution_step()

    # step by step:

    def InitializeSolutionStep(self):
        self._processes_execute_initialize_solution_step()
        self.solvers[self.solver_counter].InitializeSolutionStep()

    def SolveSolutionStep(self):
        self.solvers[self.solver_counter].SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.solvers[self.solver_counter].FinalizeSolutionStep()

        if(self.solver_counter == len(self.solvers)):
            self.solver_counter = 0
        else:
            self.solver_counter += 1
            
        self._processes_execute_finalize_solution_step()

    #### Solver internal methods ####

    #
    def _processes_execute_initialize(self):
        for process in self.list_of_processes:
            process.ExecuteInitialize()

    #
    def _processes_execute_initialize_solution_step(self):
        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

    #
    def _processes_execute_finalize_solution_step(self):
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()
