from importlib import import_module

import KratosMultiphysics as Kratos
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
            "analysis_settings": {},
            "analysis_output_settings": {
                     "nodal_solution_step_data_variables": [],
                     "nodal_data_value_variables"        : [],
                     "element_data_value_variables"      : [],
                     "condition_data_value_variables"    : []
            }
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

        self.model_parts: 'list[Kratos.ModelPart]' = []
        analysis_full_module = f"{analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(analysis_type)}"
        self.analysis: AnalysisStage = getattr(import_module(analysis_full_module), analysis_type)(self.model, analysis_settings.Clone())

        self.nodal_solution_step_data_variables = []
        if self.parameters["analysis_output_settings"].Has("nodal_solution_step_data_variables"):
            self.nodal_solution_step_data_variables = self.parameters["analysis_output_settings"]["nodal_solution_step_data_variables"]
        self.nodal_data_value_variables = []
        if self.parameters["analysis_output_settings"].Has("nodal_data_value_variables"):
            self.nodal_data_value_variables = self.parameters["analysis_output_settings"]["nodal_data_value_variables"]
        self.element_data_value_variables = []
        if self.parameters["analysis_output_settings"].Has("element_data_value_variables"):
            self.element_data_value_variables = self.parameters["analysis_output_settings"]["element_data_value_variables"]
        self.condition_data_value_variables = []
        if self.parameters["analysis_output_settings"].Has("condition_data_value_variables"):
            self.condition_data_value_variables = self.parameters["analysis_output_settings"]["condition_data_value_variables"]


    def GetAnalysisModelPart(self):
        return self.analysis._GetSolver().GetComputingModelPart()

    def Initialize(self):
        self.analysis.Initialize()

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

        optimization_problem_data = self.optimization_problem.GetProblemDataContainer()
        if optimization_problem_data.HasValue("requested_container_expression_outputs") and self in optimization_problem_data["requested_container_expression_outputs"]:
            self._OutputAnalysisData()

    def _OutputAnalysisData(self):
        unbuffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        for model_part in self.model_parts:
            for nodal_data in self.nodal_solution_step_data_variables:
                nodal_field = Kratos.Expression.NodalExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(nodal_field, Kratos.KratosGlobals.GetVariable(nodal_data.GetString()), True)
                unbuffered_data.SetValue(nodal_data.GetString(), nodal_field.Clone(), overwrite=True)
            for nodal_data in self.nodal_data_value_variables:
                nodal_field = Kratos.Expression.NodalExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(nodal_field, Kratos.KratosGlobals.GetVariable(nodal_data.GetString()), False)
                unbuffered_data.SetValue(nodal_data.GetString(), nodal_field.Clone(), overwrite=True)
            for elem_data in self.element_data_value_variables:
                elem_field = Kratos.Expression.ElementExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(elem_field, Kratos.KratosGlobals.GetVariable(elem_data.GetString()))
                unbuffered_data.SetValue(elem_data.GetString(), elem_field.Clone(), overwrite=True)
            for cond_data in self.condition_data_value_variables:
                cond_field = Kratos.Expression.ConditionExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(cond_field, Kratos.KratosGlobals.GetVariable(cond_data.GetString()))
                unbuffered_data.SetValue(cond_data.GetString(), cond_field.Clone(), overwrite=True)





