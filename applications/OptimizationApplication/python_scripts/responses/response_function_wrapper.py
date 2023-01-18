import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

class ResponseFunctionWrapper(OptimizationRoutine):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "name"     : "",
            "module"   : "KratosMultiphysics.OptimizationApplication.responses",
            "type"     : "",
            "settings" : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.__name = parameters["name"].GetString()
        if self.__name == "":
            raise RuntimeError(f"\"name\" field cannot be empty in response functions [ response function settings: {str(parameters)} ].")

        self.__response_function: ResponseFunction = Factory(parameters["module"].GetString(), parameters["type"].GetString(), model, parameters["settings"], optimization_info, ResponseFunction)

        # initialize containers to avoid multiple calculations in same iteration
        self.__response_value = None

    def Check(self):
        self.__response_function.Check()

    def Initialize(self):
        self.__response_function.Initialize()

    def InitializeSolutionStep(self):
        # initialize containers to avoid multiple calculations in same iteration
        self.__response_value = None

        self.__response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.__response_function.FinalizeSolutionStep()

    def Finalize(self):
        self.__response_function.Finalize()

    def GetName(self) -> str:
        return self.__name

    def CalculateValue(self):
        if self.__response_value is None:
            self.__response_value = self.__response_function.CalculateValue()
        return self.__response_value

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        self.__response_function.CalculateSensitivity(sensitivity_variable, sensitivity_container)

