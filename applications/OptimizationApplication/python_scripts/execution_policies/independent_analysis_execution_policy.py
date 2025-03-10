from typing import Union
from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.project import Project
from KratosMultiphysics.orchestrators.orchestrator import Orchestrator
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
            "analysis_module"         : "KratosMultiphysics",
            "analysis_type"           : "",
            "analysis_model_part_name": "",
            "analysis_settings"       : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]
        self.analysis_model_part_name = parameters["analysis_model_part_name"].GetString()

        if self.analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{self.analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"


    def Initialize(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def Execute(self):
        analysis_type = getattr(import_module(self.analysis_full_module), self.analysis_type)
        if AnalysisStage in analysis_type.mro():
            # the analysis type is derrived from AnalysisStage
            self.current_analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.model, self.analysis_settings.Clone())
        elif Orchestrator in analysis_type.mro():
            # the analysis type is derrive from the Orchestrator
            project = Project(self.analysis_settings.Clone())
            self.current_analysis: Orchestrator = getattr(import_module(self.analysis_full_module), self.analysis_type)(project)

        import pdb
        pdb.set_trace()
        self.current_analysis.Run()

    def GetAnalysisModelPart(self):
        if self.current_analysis is not None:
            if isinstance(self.current_analysis, AnalysisStage):
                if self.analysis_model_part_name == "" or self.analysis_model_part_name == self.current_analysis._GetSolver().GetComputingModelPart().FullName():
                    return self.current_analysis._GetSolver().GetComputingModelPart()
                else:
                    raise RuntimeError(f"The specified analysis model part name mismatch [ specified analysis model part name = {self.analysis_model_part_name}, used analysis model part name = {self.current_analysis._GetSolver().GetComputingModelPart().FullName()} ].")
            elif isinstance(self.current_analysis, Orchestrator):
                return self.current_analysis.GetProject().GetModel()[self.analysis_model_part_name]
            else:
                raise RuntimeError(f"Unsupported analysis type = {self.current_analysis}.")
        else:
            raise RuntimeError("The analysis is not run yet.")
