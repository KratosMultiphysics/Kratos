import KratosMultiphysics as Kratos
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import FileLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class ExecutionPolicyDecorator(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_parameters = Kratos.Parameters("""{
            "name"           : "",
            "type"           : "please_provide_python_module_name",
            "module"         : "KratosMultiphysics.OptimizationApplication.execution_policies",
            "pre_operations" : [],
            "post_operations": [],
            "log_in_file"    : false,
            "log_file_name"  : "",
            "settings"       : {}
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        super().__init__(parameters["name"].GetString())

        self.__log_file_name = parameters["log_file_name"].GetString()
        self.__log_in_file = parameters["log_in_file"].GetBool()

        if self.__log_in_file and self.__log_file_name == "":
            raise RuntimeError("Logger file name is empty for execution policy. Please provide \"log_file_name\" or make \"log_in_file\" to false.")

        factory = KratosModelParametersFactory(model)

        # create operations
        self.__list_of_pre_operations: 'list[Kratos.Operation]' = factory.ConstructListOfItems(parameters["pre_operations"])
        self.__list_of_post_operations: 'list[Kratos.Operation]' = factory.ConstructListOfItems(parameters["post_operations"])

        # create execution policy
        self.__execution_policy: ExecutionPolicy = OptimizationComponentFactory(model, parameters, optimization_problem)

    def GetExecutionPolicy(self):
        return self.__execution_policy

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.__execution_policy.GetAnalysisModelPart()

    def Initialize(self):
        self.__execution_policy.Initialize()

    def Check(self):
        self.__execution_policy.Check()

    def Finalize(self):
        self.__execution_policy.Finalize()

    def Execute(self):
        if self.__log_in_file:
            with FileLogger(self.__log_file_name):
                self.__ExecuteWithoutFileLogger()
        else:
            self.__ExecuteWithoutFileLogger()

    def __ExecuteWithoutFileLogger(self):
        CallOnAll(self.__list_of_pre_operations, Kratos.Operation.Execute)
        self.__execution_policy.Execute()
        CallOnAll(self.__list_of_post_operations, Kratos.Operation.Execute)
