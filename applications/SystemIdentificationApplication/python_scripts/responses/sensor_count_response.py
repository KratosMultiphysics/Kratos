import numpy
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorCountResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorCountResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorCountResponse(parameters["name"].GetString(), model, parameters["settings"])


class SensorCountResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name": ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorCountResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        nodal_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, KratosSI.SENSOR_STATUS)
        nodal_ta.CollectData()
        return numpy.sum(nodal_ta.data)

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        # make everything zeros
        for physical_variable, combined_ta in physical_variable_combined_tensor_adaptor.items():
            if physical_variable != KratosSI.SENSOR_STATUS:
                raise RuntimeError(f"Unsupported variable = {physical_variable.Name()}.")

            if len(combined_ta.GetTensorAdaptors()) != 1:
                raise RuntimeError(f"Currently only supports sensitivities for one model part.")

            if combined_ta.GetTensorAdaptors()[0].GetContainer() != self.model_part.Nodes:
                raise RuntimeError("Tensor adaptor container and mask status container mismatch.")

            combined_ta.GetTensorAdaptors()[0].data[:] = 1.0
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, copy=False).CollectData()

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"