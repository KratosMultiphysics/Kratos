from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import solverprint, bold

def CreateSolverWrapper(model, settings, solver_name):
    raise Exception('"KratosBaseWrapper" is a baseclass and cannot be used directly!')

class KratosBaseWrapper(CoSimulationSolverWrapper):
    def __init__(self, model, settings, solver_name):
        super(KratosBaseWrapper, self).__init__(model, settings, solver_name)

        input_file_name = self.settings["settings"]["input_file"].GetString()
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            self.project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    def Initialize(self):
        self._GetAnalysisStage().Initialize()

    def Finalize(self):
        self._GetAnalysisStage().Finalize()

    def AdvanceInTime(self, current_time):
        new_time = self._GetAnalysisStage()._GetSolver().AdvanceInTime(current_time)
        self._GetAnalysisStage().time = new_time # only needed to print the time correctly
        self.delta_time = new_time - current_time
        return new_time

    def Predict(self):
        self._GetAnalysisStage()._GetSolver().Predict()

    def InitializeSolutionStep(self):
        self._GetAnalysisStage().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetAnalysisStage().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetAnalysisStage().OutputSolutionStep()

    def SolveSolutionStep(self):
        self._GetAnalysisStage()._GetSolver().SolveSolutionStep()

    def GetDeltaTime(self):
        if not hasattr(self, 'delta_time'):
            raise Exception("DeltaTime can only be querried after it has been computed at least once")
        return self.delta_time

    def _GetAnalysisStage(self):
        if not hasattr(self, '_analysis_stage'):
            self._analysis_stage = self._CreateAnalysisStage()
        return self._analysis_stage

    def _CreateAnalysisStage(self):
        raise Exception("Creation of the AnalysisStage must be implemented in the derived class!")

    def _GetParallelType(self):
        raise Exception("Returning the type of parallelism must be implemented in the derived class!")


    def PrintInfo(self):
        solverprint("KratosSolver", bold(self._Name()))
        ## TODO print additional stuff with higher echo-level

    def IsDistributed(self):
        return (self._GetParallelType() == "MPI")
