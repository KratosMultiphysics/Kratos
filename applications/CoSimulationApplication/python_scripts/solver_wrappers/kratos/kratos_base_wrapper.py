from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_name):
    raise Exception('"KratosBaseWrapper" is a baseclass and cannot be used directly!')

class KratosBaseWrapper(CoSimulationSolverWrapper):
    """This class serves as basis for the kratos-wrappers
    It uses the AnalysisStage as interface to Kratos
    """
    def __init__(self, settings, solver_name):
        super(KratosBaseWrapper, self).__init__(settings, solver_name)

        input_file_name = self.settings["settings"]["input_file"].GetString()
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            self.project_parameters = KM.Parameters(parameter_file.read())

    def Initialize(self):
        self._GetAnalysisStage().Initialize()
        super(KratosBaseWrapper, self).Initialize()

    def Finalize(self):
        self._GetAnalysisStage().Finalize()

    def AdvanceInTime(self, current_time):
        new_time = self._GetAnalysisStage()._GetSolver().AdvanceInTime(current_time)
        self._GetAnalysisStage().time = new_time # only needed to print the time correctly
        return new_time

    def InitializeSolutionStep(self):
        self._GetAnalysisStage().InitializeSolutionStep()

    def Predict(self):
        self._GetAnalysisStage()._GetSolver().Predict()

    def SolveSolutionStep(self):
        self._GetAnalysisStage()._GetSolver().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetAnalysisStage().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetAnalysisStage().OutputSolutionStep()


    def _GetAnalysisStage(self):
        if not hasattr(self, '_analysis_stage'):
            self._analysis_stage = self._CreateAnalysisStage()
        return self._analysis_stage

    def _CreateAnalysisStage(self):
        raise Exception("Creation of the AnalysisStage must be implemented in the derived class!")


    def PrintInfo(self):
        cs_tools.cs_print_info("KratosSolver", self._Name())
        ## TODO print additional stuff with higher echo-level
