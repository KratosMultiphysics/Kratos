# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import FileLogger

def RetrieveClass(module_full_path_with_class_name: str):
    data = module_full_path_with_class_name.split(".")
    if len(data) < 2:
        raise RuntimeError(f"Please provide module name with the class name. [ provided name = \"{module_full_path_with_class_name}\" ]")

    module = import_module(".".join(data[:-1]))
    return getattr(module, data[-1])

class ExecutionPolicyWrapper:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_parameters = Kratos.Parameters("""{
            "name"                     : "",
            "execution_policy_type"    : "EXECUTION_POLICY_MODULE_FULL_NAME.EXECUTION_POLICY_CLASS_NAME",
            "execution_policy_settings": {},
            "pre_operations"           : [],
            "post_operations"          : [],
            "log_in_file"              : false,
            "log_file_name"            : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__log_file_name = parameters["log_file_name"].GetString()
        self.__log_in_file = parameters["log_in_file"].GetBool()
        self.__name = parameters["name"].GetString()

        if self.__name == "":
            raise RuntimeError(f"Execution policies should be given a non-empty name. Followings are the corresponding execution policy settings:\n{parameters}")

        if self.__log_in_file and self.__log_file_name == "":
            raise RuntimeError(f"Logger file name is empty for execution policy. Please provide \"log_file_name\" or make \"log_in_file\" to false.")

        factory = KratosModelParametersFactory(model)

        # create operations
        self.__list_of_pre_operations = factory.ConstructListOfItems(parameters["pre_operations"])
        self.__list_of_post_operations = factory.ConstructListOfItems(parameters["post_operations"])

        # create execution policy
        execution_policy_type = RetrieveClass(parameters["execution_policy_type"].GetString())
        self.__execution_policy = execution_policy_type(model, parameters["execution_policy_settings"])

    def Initialize(self, optimization_info: dict):
        self.__execution_policy.Initialize(optimization_info)

    def Execute(self, optimization_info: dict):
        if self.__log_in_file:
            with FileLogger(self.__log_file_name):
                self.__ExecuteWithoutFileLogger(optimization_info)
        else:
            self.__ExecuteWithoutFileLogger(optimization_info)

    def GetExecutionPolicy(self):
        return self.__execution_policy

    def GetExecutionPolicyName(self):
        return self.__name

    def __ExecuteWithoutFileLogger(self, optimization_info: dict):
        self.__RunPreExecutionPolicyOperations()
        self.__execution_policy.Execute(optimization_info)
        self.__RunPostExecutionPolicyOperations()

    def __RunPreExecutionPolicyOperations(self):
        for operation in self.__list_of_pre_operations:
            operation.Execute()

    def __RunPostExecutionPolicyOperations(self):
        for operation in self.__list_of_post_operations:
            operation.Execute()

