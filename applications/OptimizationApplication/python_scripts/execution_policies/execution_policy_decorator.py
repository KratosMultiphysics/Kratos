import KratosMultiphysics as Kratos
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import FileLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationProcessFactory
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class ExecutionPolicyDecorator(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationProblem):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "name"           : "",
            "module"         : "KratosMultiphysics.OptimizationApplication.execution_policies",
            "type"           : "PleaseProvideClassName",
            "settings"       : {},
            "pre_operations" : [],
            "post_operations": [],
            "log_in_file"    : false,
            "log_file_name"  : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__log_file_name = parameters["log_file_name"].GetString()
        self.__log_in_file = parameters["log_in_file"].GetBool()
        self.__name = parameters["name"].GetString()

        if self.__name == "":
            raise RuntimeError(f"Execution policies should be given a non-empty name. Followings are the corresponding execution policy settings:\n{parameters}")

        if self.__log_in_file and self.__log_file_name == "":
            raise RuntimeError("Logger file name is empty for execution policy. Please provide \"log_file_name\" or make \"log_in_file\" to false.")

        factory = KratosModelParametersFactory(model)

        # create operations
        self.__list_of_pre_operations: 'list[Kratos.Operation]' = factory.ConstructListOfItems(parameters["pre_operations"])
        self.__list_of_post_operations: 'list[Kratos.Operation]' = factory.ConstructListOfItems(parameters["post_operations"])

        # create execution policy
        self.__execution_policy: ExecutionPolicy = OptimizationProcessFactory(parameters["module"].GetString(), parameters["type"].GetString(), model, parameters["settings"], optimization_info ,ExecutionPolicy)

    def ExecuteInitialize(self):
        self.__execution_policy.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.__execution_policy.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.__execution_policy.ExecuteInitializeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.__execution_policy.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.__execution_policy.ExecuteAfterOutputStep()

    def ExecuteFinalizeSolutionStep(self):
        self.__execution_policy.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        self.__execution_policy.ExecuteFinalize()

    def Execute(self):
        if self.__log_in_file:
            with FileLogger(self.__log_file_name):
                self.__ExecuteWithoutFileLogger()
        else:
            self.__ExecuteWithoutFileLogger()

    def GetExecutionPolicy(self):
        return self.__execution_policy

    def GetExecutionPolicyName(self):
        return self.__name

    def __ExecuteWithoutFileLogger(self):
        self.__RunPreExecutionPolicyOperations()
        self.__execution_policy.Execute()
        self.__RunPostExecutionPolicyOperations()

    def __RunPreExecutionPolicyOperations(self):
        for operation in self.__list_of_pre_operations:
            operation.Execute()

    def __RunPostExecutionPolicyOperations(self):
        for operation in self.__list_of_post_operations:
            operation.Execute()

