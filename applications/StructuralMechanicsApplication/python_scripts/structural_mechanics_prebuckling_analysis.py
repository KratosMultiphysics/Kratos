from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class StructuralMechanicsPrebucklingAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters):
        solver_settings = project_parameters["solver_settings"]
        super(StructuralMechanicsPrebucklingAnalysis, self).__init__(model, project_parameters)
    
    '''Break Solution Loop when Buckling Analysis is converged.
    End Time Parameter can be considered as maximal iteration number
    '''
    def RunSolutionLoop(self):
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            if( self._GetSolver().get_mechanical_solution_strategy().GetSolutionFoundFlag() ):
                break

    def FinalizeSolutionStep(self):
        ''' This function is overriden to postprocess eigenvalues only every second load step.
        Print eigenvalues after erery small loadstep'''
        self._GetSolver().FinalizeSolutionStep()
        for process in self._GetListOfProcesses():
            if( process.__class__.__name__ != "PostprocessEigenvaluesProcess" ):
                process.ExecuteFinalizeSolutionStep() 
            elif ( (self.time % 2 == 0) & (self.time > 0 ) ):
                process.ExecuteFinalizeSolutionStep() 
                  
    def OutputSolutionStep(self):
        ''' This function is overriden to print output only every second load step.
        Print Output after every big loadstep'''
        is_output_step = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                is_output_step = True
                break

        if is_output_step: # at least one of the output processes will print output
            for process in self._GetListOfProcesses():
                process.ExecuteBeforeOutputStep()
            
            for output_process in self._GetListOfOutputProcesses():
                if( output_process.IsOutputStep() & (self.time % 2 == 1) ):
                    output_process.PrintOutput()

            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("PrebucklingAnalysis", warn_msg)
            