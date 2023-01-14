from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import RetrieveObject
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import ContainerEnum
from KratosMultiphysics.OptimizationApplication.responses.response_function import GetSensitivityContainer

class ResponseFunctionBaseWrapper(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = self.GetDefaultParameters()
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.name = parameters["name"].GetString()

        if self.name == "":
            raise RuntimeError(f"Empty \"name\" fields are not allowed in response functions. Provided response settings: \n{str(parameters)}")

        # create the response function
        response_function_settings = Kratos.Parameters("""{}""")
        response_function_settings.AddString("module", parameters["module"].GetString())
        response_function_settings.AddString("type", parameters["type"].GetString())
        response_function_settings.AddValue("settings", parameters["settings"])
        self.response_function: ResponseFunction = RetrieveObject(self.model, response_function_settings, optimization_info, ResponseFunction)

        # response storage
        self.response_value = None
        self.response_sensitivities = {}

    def Initialize(self):
        self.response_function.Initialize()

    def InitializeSolutionStep(self):
        # new design iteration, hence resetting values
        self.response_value = None
        self.response_sensitivities = {}

        self.response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.response_function.FinalizeSolutionStep()

    def Finalize(self):
        self.response_function.Finalize()

    def GetName(self):
        return self.name

    def GetValue(self):
        if self.response_value is None:
            self.response_value = self.response_function.CalculateValue()

        return self.response_value

    def GetSensitivity(self, sensitivity_variable, sensitivity_model_part: Kratos.ModelPart, sensitivity_container_type: ContainerEnum):
        if not (sensitivity_variable, sensitivity_model_part, sensitivity_container_type) in self.response_sensitivities.keys():
            # calculate the sensitivity
            self.response_function.CalculateSensitivity(sensitivity_variable, sensitivity_model_part, sensitivity_container_type)

            # transfer sensitivity values from the model part to standalone matrix/vector
            # so the same sensitivity variable can be used to store sensitivities
            # w.r.t. another response function

            container = GetSensitivityContainer(sensitivity_model_part, sensitivity_container_type)
            values = Kratos.Vector()
            domain_size = sensitivity_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
            KratosOA.OptimizationVariableUtils.GetContainerVariableToVector(container, sensitivity_variable, domain_size, values)
            self.response_sensitivities[sensitivity_variable, sensitivity_model_part, sensitivity_container_type] = values

        return self.response_sensitivities[sensitivity_variable, sensitivity_model_part, sensitivity_container_type]

    def GetStandardizedValue(self):
        return self._StandardizeValue(self.GetValue())

    def GetStandardizedSensitivity(self, sensitivity_variable, sensitivity_model_part: Kratos.ModelPart, sensitivity_container_type: ContainerEnum):
        return self._StandardizeSensitivity(self.GetSensitivity(sensitivity_variable, sensitivity_model_part, sensitivity_container_type))

    @abstractmethod
    def _StandardizeValue(self, value):
        pass

    @abstractmethod
    def _StandardizeSensitivity(self, sensitivities):
        pass

    @abstractmethod
    def GetDefaultParameters(self) -> Kratos.Parameters:
        pass

    @abstractmethod
    def IsActive(self) -> bool:
        pass

    @abstractmethod
    def GetInfo(self, spacing: str) -> str:
        pass


class ObjectiveResponseFunctionWrapper(ResponseFunctionBaseWrapper):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        self.scaling = parameters["scaling"].GetDouble()

        self.objective = parameters["objective"].GetString()
        match self.objective:
            case "minimization":
                self.standardization_value = 1.0
            case "maximization":
                self.standardization_value = -1.0
            case _:
                raise RuntimeError(f"{self.name} requesting unsupported objective {self.objective}. Supported objectives are: \n\tminimization\n\tmaximization")

    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "name"     : "",
            "module"   : "",
            "type"     : "",
            "objective": "",
            "scaling"  : 1.0,
            "settings" : {}
        }""")

    def IsActive(self) -> bool:
        return True

    def _StandardizeValue(self, value):
        return value * self.standardization_value * self.scaling

    def _StandardizeSensitivity(self, sensitivities):
        return sensitivities * self.standardization_value * self.scaling

    def GetInfo(self, spacing: str) -> str:
        info  = f"{spacing}response name              : {self.name}\n"
        info += f"{spacing}response type              : {self.response_function.__class__.__name__}\n"
        info += f"{spacing}response objective         : {self.objective}\n"
        info += f"{spacing}response value             : {self.GetValue()}\n"
        info += f"{spacing}response value             : {self.GetValue()}\n"
        info += f"{spacing}standardized response value: {self.GetStandardizedValue()}\n"
        return info


class ConstraintResponseFunctionWrapper(ResponseFunctionBaseWrapper):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        self.scaling = parameters["scaling"].GetDouble()
        self.ref_scaling = parameters["ref_scaling"].GetDouble()

        self.ref_type = parameters["ref_type"].GetString()
        match self.ref_type:
            case "initial":
                self.reference_value = None
            case "specified":
                self.reference_value = parameters["ref_value"].GetDouble()
            case _:
                raise RuntimeError(f"Provided \"reference_type\" = {self.ref_type} is not supported for response function with \"name\" = {self.name}. Followings are supported options: \n\tinitial\n\tspecified")

        self.constraint = parameters["constraint"].GetString()
        match self.constraint:
            case "<" | "<=" | "=":
                self.standardization_value = 1.0
            case ">" | ">=":
                self.standardization_value = -1.0
            case _:
                raise RuntimeError(f"Provided \"constraint\" = {self.constraint} is not supported in response function with \"name\" = {self.name}. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

    def GetDefaultParameters(self) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "name"       : "",
            "module"     : "",
            "type"       : "",
            "scaling"    : 1.0,
            "constraint" : "",
            "ref_type"   : "",
            "ref_value"  : 0.0,
            "ref_scaling": 1.0,
            "settings"   : {}
        }""")

    def IsActive(self) -> bool:
        return (self.constraint == "=") or self.GetStandardizedValue() >= 0.0

    def _StandardizeValue(self, value):
        if self.reference_value is None:
            self.reference_value = self.GetValue()
        return self.standardization_value * (value - self.reference_value * self.ref_scaling) * self.scaling

    def _StandardizeSensitivity(self, sensitivities):
        return sensitivities * self.standardization_value * self.scaling

    def GetInfo(self, spacing: str) -> str:
        info  = f"{spacing}constraint name                  : {self.name}\n"
        info += f"{spacing}constraint type                  : {self.response_function.__class__.__name__}\n"
        info += f"{spacing}constraint constraint            : {self.ref_type}\n"
        info += f"{spacing}constraint value                 : {self.GetValue()}\n"
        info += f"{spacing}constraint value                 : {self.GetValue()}\n"
        info += f"{spacing}constraint scaled reference value: {self.reference_value * self.ref_scaling}\n"
        info += f"{spacing}standardized constraint value: {self.GetStandardizedValue()}\n"
        return info


def CreateResponseFunctionWrapper(model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
    if parameters.Has("objective"):
        return ObjectiveResponseFunctionWrapper(model, parameters, optimization_info)
    elif parameters.Has("constraint"):
        return ConstraintResponseFunctionWrapper(model, parameters, optimization_info)
    else:
        raise RuntimeError(f"Either \"objective\" or \"constraint\" should be present to differentiate responses in response settings. The provided response settinsg: \n{str(parameters)}")