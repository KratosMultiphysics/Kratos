from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as tools

# Importing the Kratos Library
try:
    import KratosMultiphysics
except ModuleNotFoundError:
    raise ModuleNotFoundError(tools.bcolors.FAIL + 'KRATOS is not available! Please ensure that Kratos is available for usage!'+ tools.bcolors.ENDC)

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_co_simulation_classes.co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
import os

class KratosBaseFieldSolver(CoSimulationBaseSolver):
    def __init__(self, solver_name, cosim_solver_settings):
        super(KratosBaseFieldSolver, self).__init__(solver_name, cosim_solver_settings)

        working_directory = self.cosim_solver_settings["settings"]["working_directory"].GetString()
        input_file_name = os.path.join(working_directory, self.cosim_solver_settings["settings"]["input_file"].GetString())
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.project_parameters = KratosMultiphysics.Parameters(parameters)

    def Initialize(self):
        super(KratosBaseFieldSolver, self).Initialize()
        self._GetAnalysisStage().Initialize()


    def Finalize(self):
        super(KratosBaseFieldSolver, self).Finalize()
        self._GetAnalysisStage().Finalize()

    def AdvanceInTime(self, current_time):
        new_time = self._GetAnalysisStage()._GetSolver().AdvanceInTime(current_time)
        self._GetAnalysisStage().time = new_time # only needed to print the time correctly
        self.delta_time = new_time - current_time
        return new_time

    def Predict(self):
        super(KratosBaseFieldSolver, self).Predict()
        self._GetAnalysisStage()._GetSolver().Predict()

    def InitializeSolutionStep(self):
        super(KratosBaseFieldSolver, self).InitializeSolutionStep()
        self._GetAnalysisStage().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(KratosBaseFieldSolver, self).FinalizeSolutionStep()
        self._GetAnalysisStage().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetAnalysisStage().OutputSolutionStep()

    def SolveSolutionStep(self):
        self._GetAnalysisStage()._GetSolver().SolveSolutionStep()

    def GetBufferSize(self):
        model_part_name = self.project_parameters["solver_settings"]["model_part_name"].GetString()
        return self.model[model_part_name].GetBufferSize()

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
        solverprint(self.lvl, "KratosSolver", bold(self._Name()))
        ## TODO print additional stuff with higher echo-level

    def Check(self):
        pass

    def _GetIOName(self):
        return "kratos"