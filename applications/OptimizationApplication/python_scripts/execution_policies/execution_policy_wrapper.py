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
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import FileLogger

def RetrieveObject(model: Kratos.Model, parameters: Kratos.Parameters):
    default_settings = Kratos.Parameters("""{
        "module"  : "",
        "type"    : "",
        "settings": {}
    }""")
    parameters.ValidateAndAssignDefaults(default_settings)

    class_name = parameters["type"].GetString()
    python_file_name = ''.join(['_' + c.lower() if c.isupper() else c for c in class_name]).lstrip('_')

    module_name = parameters["module"].GetString()
    if module_name != "":
        module_name += "."
    module_name += python_file_name

    module = import_module(module_name)
    return getattr(module, class_name)(model, parameters["settings"])

class ExecutionPolicyWrapper:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_parameters = Kratos.Parameters("""{
            "name"                     : "",
            "execution_policy_settings": {
                "module"  : "KratosMultiphysics.OptimizationApplication.execution_policies",
                "type"    : "PleaseProvideClassName",
                "settings": {}
            },
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
        parameters["execution_policy_settings"].ValidateAndAssignDefaults(default_parameters["execution_policy_settings"])
        self.__execution_policy = RetrieveObject(model, parameters["execution_policy_settings"])
        if not isinstance(self.__execution_policy, ExecutionPolicy):
            raise RuntimeError(f"{self.__execution_policy.__class__.__name__} is not derrived from ExecutionPolicy.")

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

