# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import RetrieveObject

class SteppingAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "analysis_settings": {
                "module"  : "",
                "type"    : "",
                "settings": {}
            }
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_parts = []
        self.__analysis = RetrieveObject(self.model, parameters["analysis_settings"], AnalysisStage)

    def Initialize(self, _: dict):
        self.__analysis.Initialize()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def InitializeIteration(self, _: dict):
        pass

    def Execute(self, _: dict):
        time_before_analysis = []
        step_before_analysis = []
        delta_time_before_analysis = []

        for model_part in self.model_parts:
            time_before_analysis.append(model_part.ProcessInfo[Kratos.TIME])
            step_before_analysis.append(model_part.ProcessInfo[Kratos.STEP])
            delta_time_before_analysis.append(model_part.ProcessInfo[Kratos.DELTA_TIME])

        # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index] - 1)
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index] - 1)
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0)

        self.__analysis.time = self.__analysis._GetSolver().AdvanceInTime(self.__analysis.time)
        self.__analysis.InitializeSolutionStep()
        self.__analysis._GetSolver().Predict()
        self.__analysis._GetSolver().SolveSolutionStep()

        self.__analysis.FinalizeSolutionStep()
        self.__analysis.OutputSolutionStep()

        # Clear results or modifications on model parts
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_analysis[index])

    def FinalizeIteration(self, _: dict):
        pass

    def Finalize(self, _: dict):
        pass

    def GetAnalysis(self, _: dict):
        return self.__analysis

