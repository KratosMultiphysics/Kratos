import math, numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"CombinatorialResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"CombinatorialResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return CombinatorialResponse(parameters["name"].GetString(), model, parameters["settings"])


class CombinatorialResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_model_part_name": "",
            "p": 4.0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.sensor_model_part = model[parameters["sensor_model_part_name"].GetString()]
        self.p = parameters["p"].GetDouble()

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosDT.SENSOR_STATUS]

    def Initialize(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        return self.sensor_model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return None

    def CalculateValue(self) -> float:
        value = 0.0
        for node in self.sensor_model_part.Nodes:
            sensor_status = node.GetValue(KratosDT.SENSOR_STATUS)
            try:
                value -= math.log(sensor_status + 1e-6) + math.log(1 + 1e-6 - sensor_status)
            except:
                raise RuntimeError(f"sensor_status: {sensor_status}, value = {value}")
        return value

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            if physical_variable == KratosDT.SENSOR_STATUS:
                v = []
                for node in self.sensor_model_part.Nodes:
                    sensor_status = node.GetValue(KratosDT.SENSOR_STATUS)
                    v.append(-(1 - 2 * sensor_status) / ((sensor_status + 1e-6) * (1 + 1e-6 - sensor_status)))
                v = numpy.array(v)
                gradient_exp = Kratos.Expression.NodalExpression(self.sensor_model_part)
                Kratos.Expression.CArrayExpressionIO.Read(gradient_exp, v)
                collective_expression.GetContainerExpressions()[0].SetExpression(gradient_exp.GetExpression())
            else:
                raise RuntimeError("Unsupported physical variable.")

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.sensor_model_part.FullName()}]"