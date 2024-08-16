from importlib import import_module
from typing import Any

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
import KratosMultiphysics.OptimizationApplication as KratosOA

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
                     "nodal_solution_step_data_variables"  : [],
                     "nodal_data_value_variables"          : [],
                     "element_data_value_variables"        : [],
                     "element_properties_value_variables"  : [],
                     "condition_data_value_variables"      : [],
                     "condition_properties_value_variables": []
            }
        }""")
        self.model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem
        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters["analysis_output_settings"].ValidateAndAssignDefaults(default_settings["analysis_output_settings"])

        analysis_module = parameters["analysis_module"].GetString()
        analysis_type = parameters["analysis_type"].GetString()
        analysis_settings = parameters["analysis_settings"]

        if analysis_module == "KratosMultiphysics":
            analysis_full_module, analysis_type = GetClassModuleFromKratos(analysis_type)
        else:
            analysis_full_module = f"{analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(analysis_type)}"

        self.model_parts: 'list[Kratos.ModelPart]' = []
        self.analysis: AnalysisStage = getattr(import_module(analysis_full_module, package=None), analysis_type)(self.model, analysis_settings.Clone())

        analysis_output_settings = self.parameters["analysis_output_settings"]
        analysis_output_settings.ValidateAndAssignDefaults(default_settings["analysis_output_settings"])
        self.nodal_solution_step_data_variables = KratosAnalysisExecutionPolicy.__GetVariablesList(analysis_output_settings["nodal_solution_step_data_variables"].GetStringArray())
        self.nodal_data_value_variables = KratosAnalysisExecutionPolicy.__GetVariablesList(analysis_output_settings["nodal_data_value_variables"].GetStringArray())
        self.element_data_value_variables = KratosAnalysisExecutionPolicy.__GetVariablesList(analysis_output_settings["element_data_value_variables"].GetStringArray())
        self.element_properties_value_variables = KratosAnalysisExecutionPolicy.__GetVariablesList(analysis_output_settings["element_properties_value_variables"].GetStringArray())
        self.condition_data_value_variables = KratosAnalysisExecutionPolicy.__GetVariablesList(analysis_output_settings["condition_data_value_variables"].GetStringArray())
        self.condition_properties_value_variables = KratosAnalysisExecutionPolicy.__GetVariablesList(analysis_output_settings["condition_properties_value_variables"].GetStringArray())

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

        self._OutputAnalysisData()

    def _OutputAnalysisData(self):
        unbuffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        for model_part in self.model_parts:
            for variable in self.nodal_solution_step_data_variables:
                nodal_field = Kratos.Expression.NodalExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(nodal_field, variable, True)
                unbuffered_data.SetValue(variable.Name(), nodal_field.Clone(), overwrite=True)
            for variable in self.nodal_data_value_variables:
                nodal_field = Kratos.Expression.NodalExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(nodal_field, variable, False)
                unbuffered_data.SetValue(variable.Name(), nodal_field.Clone(), overwrite=True)
            for variable in self.element_data_value_variables:
                elem_field = Kratos.Expression.ElementExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(elem_field, variable)
                unbuffered_data.SetValue(variable.Name(), elem_field.Clone(), overwrite=True)
            for variable in self.element_properties_value_variables:
                elem_field = Kratos.Expression.ElementExpression(model_part)
                KratosOA.PropertiesVariableExpressionIO.Read(elem_field, variable)
                unbuffered_data.SetValue(variable.Name(), elem_field.Clone(), overwrite=True)                
            for variable in self.condition_data_value_variables:
                cond_field = Kratos.Expression.ConditionExpression(model_part)
                Kratos.Expression.VariableExpressionIO.Read(cond_field, variable)
                unbuffered_data.SetValue(variable.Name(), cond_field.Clone(), overwrite=True)
            for variable in self.condition_properties_value_variables:
                cond_field = Kratos.Expression.ConditionExpression(model_part)
                KratosOA.PropertiesVariableExpressionIO.Read(cond_field, variable)
                unbuffered_data.SetValue(variable.Name(), cond_field.Clone(), overwrite=True)                

    @staticmethod
    def __GetVariablesList(variable_names_list: 'list[str]') -> 'list[Any]':
        return [Kratos.KratosGlobals.GetVariable(variable_name) for variable_name in variable_names_list]



