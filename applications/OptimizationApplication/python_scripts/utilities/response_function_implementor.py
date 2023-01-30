from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

class ResponseFunctionImplementor(ABC):
    def __init__(self, response_name: str, optimization_info: OptimizationInfo):
        self.__response_function: ResponseFunction = optimization_info.GetOptimizationRoutine(ResponseFunction, response_name)
        self.__name = response_name

        self.__response_value = None
        self.__response_sensitivities = {}
        self.__initial_response_value = None

    def ResetResponseData(self):
        self.__response_value = None
        self.__response_sensitivities = {}

    def GetName(self) -> str:
        return self.__name

    def CalculateValue(self) -> float:
        if self.__response_value is None:
            self.__response_value = self.__response_function.CalculateValue()

        return self.__response_value

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        sensitivity_key = (sensitivity_variable, sensitivity_container.GetModelPart(), sensitivity_container.GetContainerTpe())
        if sensitivity_key not in self.__response_sensitivities.keys():
            self.__response_function.CalculateSensitivity(sensitivity_variable, sensitivity_container)
            self.__response_sensitivities[sensitivity_key] = sensitivity_container.Clone()
        else:
            sensitivity_container = self.__response_sensitivities[sensitivity_key]

    def _GetInitialResponseValue(self):
        if self.__initial_response_value is None:
            self.__initial_response_value = self.CalculateValue()

        return self.__initial_response_value

    @abstractmethod
    def CalculateStandardizedValue(self) -> float:
        pass

    @abstractmethod
    def CalculateStandardizedSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        pass

    @abstractmethod
    def GetResponseInfo(self) -> str:
        pass

class ObjectiveResponseFunctionImplementor(ResponseFunctionImplementor):
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__scaling = parameters["scaling"].GetDouble()

        self.__objective_type = parameters["type"].GetString()
        if self.__objective_type == "minimization":
            self.__scaling *= 1.0
        elif self.__objective_type == "maximization":
            self.__scaling *= -1.0
        else:
            raise RuntimeError(f"Requesting unsupported type {self.__objective_type} for objective response function. Supported types are: \n\tminimization\n\tmaximization")

        super().__init__(parameters["response_name"].GetString(), optimization_info)

    def CalculateStandardizedValue(self) -> float:
        return self.CalculateValue() * self.__scaling

    def CalculateStandardizedSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        self.CalculateSensitivity(sensitivity_variable, sensitivity_container)
        sensitivity_container *= self.__scaling

    def GetResponseInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t type          : {self.__objective_type}"
        msg += f"\n\t\t value         : {self.CalculateValue()}"
        msg += f"\n\t\t abs_change [%]: {(self.CalculateValue() / self._GetInitialResponseValue() - 1.0) * 100.0}"
        return msg

class ConstraintResponseFunctionImplementor(ResponseFunctionImplementor):
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0,
            "ref_type"     : "",
            "ref_value"    : 0.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__scaling = parameters["scaling"].GetDouble()

        self.__ref_type = parameters["ref_type"].GetString()
        if self.__ref_type == "initial_value":
            self.__reference_value = None
        elif self.__ref_type == "specified_value":
            self.__reference_value = parameters["ref_value"].GetDouble()
        else:
            raise RuntimeError(f"Provided \"reference_type\" = {self.__ref_type} is not supported for constraint response functions. Followings are supported options: \n\tinitial_value\n\tspecified_value")

        self.__constraint_type = parameters["type"].GetString()
        if self.__constraint_type in ["<", "<=", "="]:
            self.__standardization_value = 1.0
        elif self.__constraint_type in [">", ">="]:
            self.__standardization_value = -1.0
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

        super().__init__(parameters["response_name"].GetString(), optimization_info)

    def CalculateStandardizedValue(self) -> float:
        if self.__reference_value is None:
            self.__reference_value = self.CalculateValue()
        return self.__standardization_value * (self.CalculateValue() * self.__scaling - self.__reference_value)

    def CalculateStandardizedSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        self.CalculateSensitivity(sensitivity_variable, sensitivity_container)
        sensitivity_container *= (self.__standardization_value * self.__scaling)

    def IsActive(self) -> bool:
        return (self.__constraint_type == "=") or self.CalculateStandardizedValue() >= 0.0

    def GetResponseInfo(self) -> str:
        msg = "\tConstraint info:"
        msg += f"\n\t\t name     : {self.GetName()}"
        msg += f"\n\t\t value    : {self.CalculateValue()}"
        msg += f"\n\t\t type     : {self.__constraint_type}"
        msg += f"\n\t\t ref_value: {self.__reference_value}"
        msg += f"\n\t\t is_active: {self.IsActive()}"
        return msg