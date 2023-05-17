from typing import Union
from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.multistage_analysis import MultistageAnalysis
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class IndependentAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _: OptimizationProblem):
        super().__init__()

        self.model = model

        default_settings = Kratos.Parameters("""{
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]

        if self.analysis_module == "KratosMultiphysics":
            self.analysis_module = GetClassModuleFromKratos(self.analysis_type)

        self.analysis_full_module = f"{self.analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

    def ExecuteInitializeSolutionStep(self) -> None:
        self.current_analysis = None

    def Execute(self):
        self.current_analysis: Union[AnalysisStage, MultistageAnalysis] = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.model, self.analysis_settings.Clone())
        self.current_analysis.Run()

    def GetAnalysisModelPart(self):
        if self.current_analysis is not None:
            return self.current_analysis._GetSolver().GetComputingModelPart()
        else:
            raise RuntimeError("The analysis is not run for current iteration.")
