from typing import Union
from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.multistage_analysis import MultistageAnalysis
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"IndependentAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"IndependentAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return IndependentAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"])

class IndependentAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

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

    def Initialize(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def Execute(self):
        self.current_analysis: Union[AnalysisStage, MultistageAnalysis] = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.model, self.analysis_settings.Clone())
        self.current_analysis.Run()

    def GetAnalysisModelPart(self):
        if self.current_analysis is not None:
            return self.current_analysis._GetSolver().GetComputingModelPart()
        else:
            raise RuntimeError("The analysis is not run yet.")
