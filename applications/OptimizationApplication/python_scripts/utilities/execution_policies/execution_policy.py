# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

from abc import ABC
from abc import abstractmethod
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

class ExecutionPolicy(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_parameters = Kratos.Parameters("""{
            "pre_operations" : [],
            "post_operations": [],
            "settings"       : {},
            "log_in_file"    : true,
            "log_file_name"  : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.parameters = parameters
        self.log_file_name = self.parameters["log_file_name"].GetString()
        self.log_in_file = self.parameters["log_in_file"].GetBool()

        if self.log_in_file and self.log_file_name == "":
            raise RuntimeError(f"Logger file name is empty for execution policy. Please provide \"log_file_name\" or make \"log_in_file\" to false.")

        factory = KratosModelParametersFactory(self.model)

        # create operations
        self.list_of_pre_operations = factory.ConstructListOfItems(self.parameters["pre_operations"])
        self.list_of_post_operations = factory.ConstructListOfItems(self.parameters["post_operations"])

    def RunPreAnalysisOperations(self):
        for operation in self.list_of_pre_operations:
            operation.Execute()

    def RunPostAnalysisOperations(self):
        for operation in self.list_of_post_operations:
            operation.Execute()

    @abstractmethod
    def RunAnalysis(self, optimization_info: dict):
        pass

    @abstractmethod
    def Initialize(self, optimization_info: dict):
        pass

    def Run(self, optimization_info: dict):
        if self.log_in_file:
            with FileLogger(self.log_file_name):
                self.__RunWithoutFileLogger(optimization_info)
        else:
            self.__RunWithoutFileLogger(optimization_info)

    def __RunWithoutFileLogger(self, optimization_info: dict):
        self.RunPreAnalysisOperations()
        self.RunAnalysis(optimization_info)
        self.RunPostAnalysisOperations()

