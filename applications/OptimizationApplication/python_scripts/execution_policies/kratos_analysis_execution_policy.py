from importlib import import_module

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"KratosAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"KratosAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return KratosAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class KratosAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "analysis_module"  : "KratosMultiphysics",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")
        self.model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem
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

    def Initialize(self):
        self.analysis.Initialize()
        exe_pol = self.optimization_problem.GetExecutionPolicy(self.GetName())
        self.un_buffered_data = ComponentDataView(exe_pol, self.optimization_problem).GetUnBufferedData()

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.analysis.Finalize()

    def Execute(self):
        self.analysis.time = self.analysis.project_parameters["problem_data"]["start_time"].GetDouble()
        for model_part in self.model_parts:
            model_part.ProcessInfo.SetValue(Kratos.STEP, 0)
            model_part.ProcessInfo.SetValue(Kratos.TIME, 0)
            model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0)
        self.analysis.RunSolutionLoop()

        for model_part in self.model_parts:
            self.displ_field = Kratos.Expression.NodalExpression(model_part)
            Kratos.Expression.VariableExpressionIO.Read(self.displ_field, Kratos.DISPLACEMENT, True)
            self.un_buffered_data.SetValue("DISPLACEMENT", self.displ_field.Clone(), overwrite=True)



