import abc, numpy
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
        raise RuntimeError(f"SensorIsolationResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorIsolationResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorIsolationResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class IsolationMethod(abc.ABC):
    @abc.abstractmethod
    def GetData(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        pass

class SensorPositionBasedIsolationMethod(IsolationMethod):
    def __init__(self, model: Kratos.Model, sensor_group_name: str, parameters: Kratos.Parameters, _):
        default_settings = Kratos.Parameters("""{
            "isolation_method": "position_based"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.sensor_group_name = sensor_group_name

    def GetData(self):
        nodal_positions = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model[self.sensor_group_name].Nodes, Kratos.Configuration.Current)
        nodal_positions.CollectData()
        return nodal_positions

class SensorMaskBasedIsolationMethod(IsolationMethod):
    def __init__(self, model: Kratos.Model, sensor_group_name: str, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "isolation_method": "mask_based",
            "sensor_mask_name": ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.optimization_problem = optimization_problem
        self.sensor_group_name = sensor_group_name
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()

    def GetData(self):
        sensor_mask_status = None
        for mask_status_controller in GetMaskStatusControllers(ComponentDataView(self.sensor_group_name, self.optimization_problem), self.sensor_mask_name):
            if isinstance(mask_status_controller, KratosSI.SensorMaskStatus):
                sensor_mask_status = mask_status_controller

        if sensor_mask_status is None:
            raise RuntimeError(f"No sensor mask status controller is found for sensor mask \"{self.sensor_mask_name}\" in sensor group \"{self.sensor_group_name}\".")

        data = numpy.array(sensor_mask_status.GetMasks().transpose(), dtype=numpy.float64)
        return Kratos.TensorAdaptors.DoubleTensorAdaptor(sensor_mask_status.GetSensorModelPart().Nodes, Kratos.DoubleNDData(data))

class SensorIsolationResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name"        : "",
            "isolation_radius"         : 0.0,
            "isolation_method_settings": {
                "isolation_method" : "position_based"
            }
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.distance_matrix: Optional[KratosSI.DistanceMatrix] = None
        self.optimization_problem = optimization_problem

        self.isolation_radius = parameters["isolation_radius"].GetDouble()
        if self.isolation_radius <= 0.0:
            raise RuntimeError(f"The isolation radius should be positive value [ isolation_radius = {self.isolation_radius} ].")

        isolation_methods_dict = {
            "position_based": SensorPositionBasedIsolationMethod,
            "mask_based"    : SensorMaskBasedIsolationMethod
        }
        isolation_method = parameters["isolation_method_settings"]["isolation_method"].GetString()
        if isolation_method in isolation_methods_dict.keys():
            self.isolation_method: IsolationMethod = isolation_methods_dict[isolation_method](self.model, self.sensor_group_name, parameters["isolation_method_settings"], self.optimization_problem)
        else:
            raise RuntimeError(f"Unsupported isolation method = \"{isolation_method}\". Followings are supported:\n\t" + "\n\t".join(isolation_methods_dict.keys()))

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        if not sensor_group_data.GetUnBufferedData().HasValue("distance_matrix"):
            self.distance_matrix = KratosSI.DistanceMatrix()
            self.distance_matrix.Update(self.isolation_method.GetData())
            sensor_group_data.GetUnBufferedData().SetValue("distance_matrix", self.distance_matrix)
        else:
            self.distance_matrix = sensor_group_data.GetUnBufferedData().GetValue("distance_matrix")

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorIsolationResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        return KratosSI.Responses.SensorIsolationResponseUtils.CalculateValue(self.model_part, self.isolation_radius, self.distance_matrix)

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        # make everything zeros
        for physical_variable, combined_ta in physical_variable_combined_tensor_adaptor.items():
            if physical_variable != KratosSI.SENSOR_STATUS:
                raise RuntimeError(f"Unsupported variable = {physical_variable.Name()}.")

            if len(combined_ta.GetTensorAdaptors()) != 1:
                raise RuntimeError(f"Currently only supports sensitivities for one model part.")

            if combined_ta.GetTensorAdaptors()[0].GetContainer() != self.model_part.Nodes:
                raise RuntimeError("Tensor adaptor container and mask status container mismatch.")

            combined_ta.GetTensorAdaptors()[0].data[:] = KratosSI.Responses.SensorIsolationResponseUtils.CalculateGradient(self.model_part, self.isolation_radius, self.distance_matrix).data
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, copy=False).CollectData()

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"