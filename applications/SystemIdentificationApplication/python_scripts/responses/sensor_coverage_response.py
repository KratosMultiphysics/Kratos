from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers

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
            "sensor_group_name": "",
            "sensor_mask_name" : "",
            "echo_level"       : 0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.echo_level = parameters["echo_level"].GetInt()

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.mask_status: 'Optional[KratosSI.SensorMaskStatus]' = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for controller in GetMaskStatusControllers(sensor_group_data, self.sensor_mask_name):
            if isinstance(controller, KratosSI.SensorMaskStatus):
                self.mask_status = controller

        if self.mask_status == None:
            raise RuntimeError(f"SensorMaskStatus controller not found for the sensor_mask_name = \"{self.sensor_mask_name}\" in the sensor_group = \"{self.sensor_group_name}\".")

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorCoverageResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        return KratosSI.Responses.SensorCoverageResponseUtils.CalculateValue(self.mask_status)

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        # make everything zeros
        for physical_variable, combined_ta in physical_variable_combined_tensor_adaptor.items():
            if physical_variable != KratosSI.SENSOR_STATUS:
                raise RuntimeError(f"Unsupported variable = {physical_variable.Name()}.")

            if len(combined_ta.GetTensorAdaptors()) != 1:
                raise RuntimeError(f"Currently only supports sensitivities for one model part.")

            if combined_ta.GetTensorAdaptors()[0].GetContainer() != self.mask_status.GetSensorModelPart().Nodes:
                raise RuntimeError("Tensor adaptor container and mask status container mismatch.")

            combined_ta.GetTensorAdaptors()[0].data[:] = KratosSI.Responses.SensorCoverageResponseUtils.CalculateGradient(self.mask_status).data
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, copy=False).CollectData()

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"