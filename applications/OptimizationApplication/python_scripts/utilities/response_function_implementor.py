from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion

class ResponseFunctionImplementor(ABC):
    def __init__(self, response_name: str, optimization_info: OptimizationInfo):
        self.__response_function: ResponseFunction = optimization_info.GetOptimizationProcess(ResponseFunction, response_name)
        self.__name = response_name

        self.__initial_response_value = None
        self._optimization_info = optimization_info
        self._key_prefix = f"problem_data/response_data/{self.__response_function.GetModelPart().FullName()}/{self.GetName()}"

    def GetName(self) -> str:
        return self.__name

    def GetValue(self, solution_step_index: int = 0) -> float:
        key = f"{self._key_prefix}/value"
        if solution_step_index == 0 and not self._optimization_info.HasValue(key):
            self._optimization_info.SetValue(key, self.__response_function.CalculateValue())
        return self._optimization_info.GetValue(key, solution_step_index)

    def GetSensitivity(self, sensitivity_variable, sensitivity_data_container: ContainerVariableDataHolderUnion, solution_step_index: int = 0):
        key = f"problem_data/response_data/{sensitivity_data_container.GetModelPart().FullName()}/{self.GetName()}/sensitivities/{sensitivity_variable.Name()}/raw"
        if solution_step_index == 0 and not self._optimization_info.HasValue(key):
            # calculate the sensitivities
            self.__response_function.CalculateSensitivity(sensitivity_variable, sensitivity_data_container.GetModelPart())

            # transfer the model part variable data to ContainerVariableDataHolder object
            sensitivity_data_container.ReadDataFromContainerVariable(sensitivity_variable)

            # store it for future use if required.
            self._optimization_info.SetValue(key, sensitivity_data_container.Clone())
        else:
            sensitivity_data_container.CopyDataFrom(self._optimization_info.GetValue(key, solution_step_index))

    def _GetInitialResponseValue(self):
        if self.__initial_response_value is None:
            self.__initial_response_value = self.GetValue()

        return self.__initial_response_value

    @abstractmethod
    def GetStandardizedValue(self, solution_step_index: int = 0) -> float:
        pass

    @abstractmethod
    def GetStandardizedSensitivity(self, sensitivity_variable, sensitivity_data_container: ContainerVariableDataHolderUnion, solution_step_index: int = 0):
        pass

    @abstractmethod
    def GetResponseType(self) -> str:
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

    def GetStandardizedValue(self, solution_step_index: int = 0) -> float:
        return self.GetValue(solution_step_index) * self.__scaling

    def GetStandardizedSensitivity(self, sensitivity_variable, sensitivity_data_container: ContainerVariableDataHolderUnion, solution_step_index: int = 0):
        self.GetSensitivity(sensitivity_variable, sensitivity_data_container, solution_step_index)
        sensitivity_data_container *= self.__scaling

    def GetResponseType(self) -> str:
        return self.__objective_type

    def GetResponseInfo(self) -> str:

        if self._optimization_info["step"] > 1:
            relative_change = (self.GetValue() / self.GetValue(1) - 1.0)
        else:
            relative_change = 0.0

        absolute_change = (self.GetValue() / self._GetInitialResponseValue() - 1.0)

        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t type          : {self.__objective_type}"
        msg += f"\n\t\t value         : {self.GetValue()}"
        msg += f"\n\t\t rel_change [%]: {relative_change * 100.0}"
        msg += f"\n\t\t abs_change [%]: {absolute_change * 100.0}"

        self._optimization_info.SetValue(f"{self._key_prefix}/relative_change", relative_change)
        self._optimization_info.SetValue(f"{self._key_prefix}/absolute_change", absolute_change)
        self._optimization_info.SetValue(f"{self._key_prefix}/type", self.__objective_type)
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

    def GetReferenceValue(self):
        if self.__reference_value is None:
            self.__reference_value = self.GetValue()

        return self.__reference_value

    def GetStandardizedValue(self, solution_step_index: int = 0) -> float:
        return self.__standardization_value * (self.GetValue(solution_step_index) * self.__scaling - self.GetReferenceValue())

    def GetStandardizedSensitivity(self, sensitivity_variable, sensitivity_data_container: ContainerVariableDataHolderUnion, solution_step_index: int = 0):
        self.GetSensitivity(sensitivity_variable, sensitivity_data_container, solution_step_index)
        sensitivity_data_container *= (self.__standardization_value * self.__scaling)

    def IsActive(self, solution_step_index: int = 0) -> bool:
        return (self.__constraint_type == "=") or self.GetStandardizedValue(solution_step_index) >= 0.0

    def GetResponseType(self) -> str:
        return self.__constraint_type

    def GetResponseInfo(self) -> str:
        if self._optimization_info["step"] > 1:
            relative_change = (self.GetValue() / self.GetValue(1) - 1.0)
        else:
            relative_change = 0.0
        absolute_change = (self.GetValue() / self._GetInitialResponseValue() - 1.0)

        if self.IsActive():
            violation = self.GetStandardizedValue() / self.GetReferenceValue()
        else:
            violation = 0.0

        msg = "\tConstraint info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t value         : {self.GetValue()}"
        msg += f"\n\t\t type          : {self.__constraint_type}"
        msg += f"\n\t\t ref_value     : {self.GetReferenceValue()}"
        msg += f"\n\t\t is_active     : {self.IsActive()}"
        msg += f"\n\t\t rel_change [%]: {relative_change * 100.0}"
        msg += f"\n\t\t abs_change [%]: {absolute_change * 100.0}"
        msg += f"\n\t\t violation  [%]: {violation * 100.0}"

        self._optimization_info.SetValue(f"{self._key_prefix}/relative_change", relative_change)
        self._optimization_info.SetValue(f"{self._key_prefix}/absolute_change", absolute_change)
        self._optimization_info.SetValue(f"{self._key_prefix}/reference_value", self.GetReferenceValue())
        self._optimization_info.SetValue(f"{self._key_prefix}/is_active", self.IsActive())
        self._optimization_info.SetValue(f"{self._key_prefix}/type", self.__constraint_type)
        self._optimization_info.SetValue(f"{self._key_prefix}/violation_ratio", violation)

        return msg