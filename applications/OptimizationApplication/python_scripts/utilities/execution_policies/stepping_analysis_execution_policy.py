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
from KratosMultiphysics.OptimizationApplication.utilities.execution_policies.execution_policy_wrapper import RetrieveClass

class SteppingAnalysisExecutionPolicy:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):

        default_settings = Kratos.Parameters("""{
            "analysis_type"      : "PLEASE_PROVIDE_FULL_MODULE_PATH.CLASS_NAME_WITH_RUN_METHOD",
            "model_part_names"   : [],
            "analysis_parameters": {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters = parameters
        self.model = model

        self.model_parts = []

        analysis_class = RetrieveClass(parameters["analysis_type"].GetString())
        self.analysis = analysis_class(self.model, parameters["analysis_parameters"])
        if not isinstance(self.analysis, AnalysisStage):
            raise RuntimeError(f"The analysis class {self.analysis.__class__.__name__} is not derrived from AnalysisStage.")

    def Initialize(self, _: dict):
        self.analysis.Initialize()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

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

        self.analysis.time = self.analysis._GetSolver().AdvanceInTime(self.analysis.time)
        self.analysis.InitializeSolutionStep()
        self.analysis._GetSolver().Predict()
        self.analysis._GetSolver().SolveSolutionStep()

        self.analysis.FinalizeSolutionStep()
        self.analysis.OutputSolutionStep()

        # Clear results or modifications on model parts
        for index, model_part in enumerate(self.model_parts):
            model_part.ProcessInfo.SetValue(Kratos.STEP, step_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.TIME, time_before_analysis[index])
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, delta_time_before_analysis[index])



