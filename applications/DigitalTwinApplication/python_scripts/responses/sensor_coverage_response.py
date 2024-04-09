from typing import Optional
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorCoverageResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorCoverageResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorCoverageResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class SensorCoverageResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "mask_model_part_name": "",
            "mask_expression_name": "YOUNG_MODULUS_SENSITIVITY_mask"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for SensorCoverageResponse. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.mask_expression_name = parameters["mask_expression_name"].GetString()
        self.mask_model_part_name = parameters["mask_model_part_name"].GetString()

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosDT.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.mask_model_part = self.model[self.mask_model_part_name]

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorCoverageResponse::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return None

    def CalculateValue(self) -> float:
        # now compute the smooth clamped values
        smooth_clamper = KratosDT.ElementSmoothClamper(0, 1)
        clamped_theta_tilde_exp = smooth_clamper.Clamp(self.__GetEntitySummation())

        domain_size_exp = Kratos.Expression.ElementExpression(self.mask_model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)

        return Kratos.Expression.Utils.InnerProduct(domain_size_exp, clamped_theta_tilde_exp) / Kratos.Expression.Utils.Sum(domain_size_exp)

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        # now compute the smooth clamped values
        smooth_clamper = KratosDT.ElementSmoothClamper(0, 1)
        clamped_theta_tilde_derivative_exp = smooth_clamper.ClampDerivative(self.__GetEntitySummation())

        domain_size_exp = Kratos.Expression.ElementExpression(self.mask_model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)

        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = ComponentDataView("sensor", self.optimization_problem).GetUnBufferedData().GetValue("list_of_sensors")
        gradient = np.array([0.0] * len(list_of_sensors))
        for i, sensor in enumerate(list_of_sensors):
            mask = sensor.GetElementExpression(self.mask_expression_name)
            gradient[i] += Kratos.Expression.Utils.InnerProduct(mask, domain_size_exp * clamped_theta_tilde_derivative_exp)

        total_domain_size = Kratos.Expression.Utils.Sum(domain_size_exp)
        Kratos.Expression.CArrayExpressionIO.Read(physical_variable_collective_expressions[KratosDT.SENSOR_STATUS].GetContainerExpressions()[0], gradient / total_domain_size)

    def __GetEntitySummation(self) -> Kratos.Expression.ElementExpression:
        n = self.mask_model_part.NumberOfElements()

        sensor_status_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(sensor_status_exp, KratosDT.SENSOR_STATUS, False)
        sensor_status_exp = sensor_status_exp.Evaluate()

        # first we calculate \tilde{\theta_j} for each element
        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = ComponentDataView("sensor", self.optimization_problem).GetUnBufferedData().GetValue("list_of_sensors")
        theta_tilde = np.array([0.0] * n)
        for i, sensor in enumerate(list_of_sensors):
            mask = sensor.GetElementExpression(self.mask_expression_name).Evaluate()
            for j in range(n):
                theta_tilde[j] += mask[j] * sensor_status_exp[i]
        theta_tilde_exp = Kratos.Expression.ElementExpression(self.mask_model_part)
        Kratos.Expression.CArrayExpressionIO.Read(theta_tilde_exp, theta_tilde)
        return theta_tilde_exp

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"