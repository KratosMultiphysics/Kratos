import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import ResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine

class ObjectiveResponseFunctionWrapper:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.optimization_info = optimization_info
        self.parameters = parameters
        self.scaling = parameters["scaling"].GetDouble()
        self.__response_function_wrapper = None

        self.objective_type = parameters["type"].GetString()
        if self.objective_type == "minimization":
            self.standardization_value = 1.0
        elif self.objective_type == "maximization":
            self.standardization_value = -1.0
        else:
            raise RuntimeError(f"Requesting unsupported type {self.objective_type} for objective response function. Supported types are: \n\tminimization\n\tmaximization")

    def GetName(self) -> str:
        return self.__GetResponseFunction().GetName()

    def CalculateValue(self) -> float:
        return self.__GetResponseFunction().CalculateValue()

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        self.__GetResponseFunction().CalculateSensitivity(sensitivity_variable, sensitivity_container)

    def CalculateStandardizedValue(self) -> float:
        return self.CalculateValue() * self.standardization_value * self.scaling

    def CalculateStandardizedSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        self.CalculateSensitivity(sensitivity_variable, sensitivity_container)
        sensitivity_container *= self.standardization_value * self.scaling

    def GetResponseInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name : {self.GetName()}"
        msg += f"\n\t\t type : {self.objective_type}"
        msg += f"\n\t\t value: {self.CalculateValue()}"
        return msg

    def __GetResponseFunction(self):
        if self.__response_function_wrapper is None:
            self.__response_function_wrapper: ResponseFunctionWrapper = self.optimization_info.GetOptimizationRoutine("ResponseFunctionWrapper", self.parameters["response_name"].GetString())
        return self.__response_function_wrapper

class ConstraintResponseFunctionWrapper:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0,
            "ref_type"     : "",
            "ref_value"    : 0.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.optimization_info = optimization_info
        self.parameters = parameters
        self.scaling = parameters["scaling"].GetDouble()
        self.__response_function_wrapper = None

        self.ref_type = parameters["ref_type"].GetString()
        if self.ref_type == "initial_value":
            self.reference_value = None
        elif self.ref_type == "specified_value":
            self.reference_value = parameters["ref_value"].GetDouble()
        else:
            raise RuntimeError(f"Provided \"reference_type\" = {self.ref_type} is not supported for constraint response functions. Followings are supported options: \n\tinitial\n\tspecified")

        self.constraint_type = parameters["type"].GetString()
        if self.constraint_type in ["<", "<=", "="]:
            self.standardization_value = 1.0
        elif self.constraint_type in [">", ">="]:
            self.standardization_value = -1.0
        else:
            raise RuntimeError(f"Provided \"type\" = {self.constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

    def GetName(self) -> str:
        return self.__GetResponseFunction().GetName()

    def CalculateValue(self) -> float:
        return self.__GetResponseFunction().CalculateValue()

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_contaienr: ContainerData):
        self.__GetResponseFunction().CalculateSensitivity(sensitivity_variable, sensitivity_contaienr)

    def CalculateStandardizedValue(self) -> float:
        if self.reference_value is None:
            self.reference_value = self.CalculateValue()
        return self.standardization_value * (self.CalculateValue() * self.scaling - self.reference_value)

    def CalculateStandardizedSensitivity(self, sensitivity_variable, sensitivity_contaienr: ContainerData):
        self.CalculateSensitivity(sensitivity_variable, sensitivity_contaienr)
        sensitivity_contaienr *= self.standardization_value * self.scaling

    def IsActive(self) -> bool:
        return (self.constraint_type == "=") or self.CalculateStandardizedValue() >= 0.0

    def GetResponseInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name  : {self.GetName()}"
        msg += f"\n\t\t type  : {self.constraint_type} {self.reference_value}"
        msg += f"\n\t\t value : {self.CalculateValue()}"
        msg += f"\n\t\t active: {self.IsActive()}"
        return msg

    def __GetResponseFunction(self):
        if self.__response_function_wrapper is None:
            self.__response_function_wrapper: ResponseFunctionWrapper = self.optimization_info.GetOptimizationRoutine("ResponseFunctionWrapper", self.parameters["response_name"].GetString())
        return self.__response_function_wrapper

class AlgorithmWrapper(OptimizationRoutine):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "name"         : "",
            "module"       : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"         : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "control_names": [],
            "objectives"   : [],
            "constraints"  : [],
            "settings"     : {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.optimization_info = optimization_info
        self.model = model
        self.parameters = parameters

        self.__name = parameters["name"].GetString()
        if self.__name == "":
            raise RuntimeError(f"An empty algorithm name is not allowed. [ Algorithm settings: {str(parameters)} ].")

        self.__algorithm: Algorithm = Factory(parameters["module"].GetString(), parameters["type"].GetString(), model, parameters["settings"], optimization_info, Algorithm)
        self.__algorithm.SetName(self.__name)

    def AddVariables(self):
        self.__algorithm.AddVariables()

    def AddDofs(self):
        self.__algorithm.AddDofs()

    def Initialize(self):
        self.__algorithm.SetObjectives([ObjectiveResponseFunctionWrapper(self.model, objective_parameters, self.optimization_info) for objective_parameters in self.parameters["objectives"]])
        self.__algorithm.SetConstraints([ConstraintResponseFunctionWrapper(self.model, objective_parameters, self.optimization_info) for objective_parameters in self.parameters["constraints"]])
        self.__algorithm.SetControllers([self.optimization_info.GetOptimizationRoutine("ControlWrapper", control_name) for control_name in self.parameters["control_names"].GetStringArray()])

        self.__algorithm.Initialize()

    def InitializeSolutionStep(self):
        self.__algorithm.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.__algorithm.FinalizeSolutionStep()

    def Finalize(self):
        self.__algorithm.Finalize()

    def GetMinimumBufferSize(self) -> int:
        return self.__algorithm.GetMinimumBufferSize()

    def SolveSolutionStep(self):
        self.__algorithm.SolveSolutionStep()

    def IsConverged(self) -> bool:
        return self.__algorithm.IsConverged()

    def GetName(self) -> str:
        return self.__name