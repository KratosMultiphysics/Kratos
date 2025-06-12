from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"SteppingAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SteppingAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return SteppingAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class SteppingAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

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
            analysis_full_module, analysis_type = GetClassModuleFromKratos(analysis_type)
        else:
            analysis_full_module = f"{analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(analysis_type)}"

        self.model_parts = []
        self.analysis: AnalysisStage = getattr(import_module(analysis_full_module), analysis_type)(self.model, analysis_settings.Clone())

    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def Initialize(self):
        self.analysis.Initialize()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

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



