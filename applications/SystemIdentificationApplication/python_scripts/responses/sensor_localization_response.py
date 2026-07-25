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
        raise RuntimeError(f"SensorLocalizationResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorLocalizationResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorLocalizationResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class DissimilarityMethod(abc.ABC):
    pass

    @abc.abstractmethod
    def ComputeDissimilarity(self) -> float:
        pass


class ConstantDissimilarityMethod(DissimilarityMethod):
    def __init__(self, _, parameters: Kratos.Parameters, __):
        default_settings = Kratos.Parameters("""{
            "dissimilarity_type" : "constant",
            "dissimilarity_value": 1.0
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)
        self.dissimilarity_value = parameters["dissimilarity_value"].GetDouble()

    def ComputeDissimilarity(self) -> float:
        return self.dissimilarity_value

class ExponentiallyDecayingDissimilarityMethod(DissimilarityMethod):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "dissimilarity_type": "exponentially_decaying",
            "initial_value"     : 1e+6,
            "final_value"       : 1e-3,
            "decay_period"      : 100
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.optimization_problem = optimization_problem

        self.initial_value = parameters["initial_value"].GetDouble()
        self.final_value = parameters["final_value"].GetDouble()
        self.beta = (numpy.log(self.initial_value) - numpy.log(self.final_value)) / parameters["decay_period"].GetInt()

    def ComputeDissimilarity(self) -> float:
        dissimilarity_value = max(self.final_value, self.initial_value * numpy.exp(-self.beta * self.optimization_problem.GetStep()))
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Using a dissimilarity value of {dissimilarity_value}.")
        return dissimilarity_value

class StatusBasedExponentiallyDecayingDissimilarityMethod(DissimilarityMethod):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "dissimilarity_type": "status_based_exponentially_decaying",
            "sensor_group_name" : "",
            "max_value"         : 1e+6,
            "min_value"         : 1e-3
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.optimization_problem = optimization_problem

        self.max_value = parameters["max_value"].GetDouble()
        self.min_value = parameters["min_value"].GetDouble()
        self.sensor_group_name = parameters["sensor_group_name"].GetString()

    def ComputeDissimilarity(self) -> float:
        sensor_status = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model[self.sensor_group_name].Nodes, KratosSI.SENSOR_STATUS)
        sensor_status.CollectData()

        max_gap = numpy.max(numpy.abs(sensor_status.data - 0.0) * numpy.abs(1.0 - sensor_status.data))

        # the max value which the max gap can take is 0.25, where all the sensor statuses are 0.5
        # the min value which the max gap can take is 0.0 where all the sensors are either 0.0 or 1.0

        dissimilarity_value = max(self.min_value, self.max_value * numpy.exp(-(numpy.log(self.min_value/self.max_value)*(max_gap-0.25)/0.25)))
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Using a dissimilarity value of {dissimilarity_value} for the sensor status maximum gap of {max_gap}.")
        return dissimilarity_value

class SensorLocalizationResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name"     : "",
            "sensor_mask_name"      : "",
            "echo_level"            : 0,
            "boltzmann_beta"        : 1.0,
            "epsilon"               : 1e-3,
            "dissimilarity_settings": {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.echo_level = parameters["echo_level"].GetInt()
        self.boltzmann_beta = parameters["boltzmann_beta"].GetDouble()
        self.epsilon = parameters["epsilon"].GetDouble()

        dissimilarity_method_dict = {
            "constant": ConstantDissimilarityMethod,
            "exponentially_decaying": ExponentiallyDecayingDissimilarityMethod,
            "status_based_exponentially_decaying": StatusBasedExponentiallyDecayingDissimilarityMethod
        }
        dissimilarity_type = parameters["dissimilarity_settings"]["dissimilarity_type"].GetString()
        if dissimilarity_type in dissimilarity_method_dict.keys():
            self.dissimilarity_method: DissimilarityMethod = dissimilarity_method_dict[dissimilarity_type](model, parameters["dissimilarity_settings"], optimization_problem)
        else:
            raise RuntimeError(f"Unsupported dissimilarity_type = \"{dissimilarity_type}\" provided. Followings are supported:\n\t" + "\n\t".join(list(dissimilarity_method_dict.keys())))

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.mask_status_kd_tree: 'Optional[KratosSI.SensorMaskStatusKDTree]' = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for controller in GetMaskStatusControllers(sensor_group_data, self.sensor_mask_name):
            if isinstance(controller, KratosSI.SensorMaskStatusKDTree):
                self.mask_status_kd_tree = controller

        if self.mask_status_kd_tree == None:
            raise RuntimeError(f"SensorMaskStatusKDTree controller not found for the sensor_mask_name = \"{self.sensor_mask_name}\" in the sensor_group = \"{self.sensor_group_name}\".")

        self.utils = KratosSI.Responses.SensorLocalizationResponseUtils(self.mask_status_kd_tree, self.boltzmann_beta, self.epsilon)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorLocalizationResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        return self.utils.CalculateValue(self.dissimilarity_method.ComputeDissimilarity())

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        # make everything zeros
        for physical_variable, combined_ta in physical_variable_combined_tensor_adaptor.items():
            if physical_variable != KratosSI.SENSOR_STATUS:
                raise RuntimeError(f"Unsupported variable = {physical_variable.Name()}.")

            if len(combined_ta.GetTensorAdaptors()) != 1:
                raise RuntimeError(f"Currently only supports sensitivities for one model part.")

            if combined_ta.GetTensorAdaptors()[0].GetContainer() != self.mask_status_kd_tree.GetSensorMaskStatus().GetSensorModelPart().Nodes:
                raise RuntimeError("Tensor adaptor container and mask status container mismatch.")

            combined_ta.GetTensorAdaptors()[0].data[:] = self.utils.CalculateGradient(self.dissimilarity_method.ComputeDissimilarity()).data
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(combined_ta, perform_collect_data_recursively=False, copy=False).CollectData()

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"