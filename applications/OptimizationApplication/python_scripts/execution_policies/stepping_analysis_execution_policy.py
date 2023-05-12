from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class SteppingAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _: OptimizationProblem):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)

        analysis_module = parameters["analysis_module"].GetString()
        analysis_type = parameters["analysis_type"].GetString()
        analysis_settings = parameters["analysis_settings"]

        if analysis_module == "KratosMultiphysics":
            analysis_module = GetClassModuleFromKratos(analysis_type)

        self.model_parts = []
        analysis_full_module = f"{analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(analysis_type)}"
        self.analysis: AnalysisStage = getattr(import_module(analysis_full_module), analysis_type)(self.model, analysis_settings.Clone())

    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def ExecuteInitialize(self):
        self.analysis.Initialize()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Execute(self):
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



