from importlib import import_module
from typing import Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"KratosCoSimAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"KratosCoSimAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return KratosCoSimAnalysisExecutionPolicy(parameters["name"].GetString(), models, parameters["settings"], optimization_problem)

class KratosCoSimAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_names" : [],
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")
        self.models = models
        self.parameters = parameters
        self.optimization_problem = optimization_problem
        self.parameters.ValidateAndAssignDefaults(default_settings)

        analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_settings = parameters["analysis_settings"]

        if analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

        self.analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(self.analysis_settings.Clone(), self.models)

        self.variable_utils = Kratos.VariableUtils()

    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def Initialize(self):
        self.analysis.Initialize()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.analysis.Finalize()

    def Execute(self):
        for k,v in self.analysis.models.items():
            model_part_names = v.GetModelPartNames()
            for model_part in model_part_names:
                mdp = v.GetModelPart(model_part)
                mdp.ProcessInfo.SetValue(Kratos.IS_RESTARTED, False)
                mdp.ProcessInfo.SetValue(Kratos.STEP, 0)
                mdp.ProcessInfo.SetValue(Kratos.TIME, 0.0)
                mdp.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.0)
        self.analysis.time = 0.0
        self.analysis.step = 0
        self.analysis.RunSolutionLoop()

    @staticmethod
    def __GetVariablesList(variable_names_list: 'list[str]') -> 'list[Any]':
        return [Kratos.KratosGlobals.GetVariable(variable_name) for variable_name in variable_names_list]



