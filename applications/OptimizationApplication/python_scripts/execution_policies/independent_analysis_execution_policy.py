# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import RetrieveObject

class IndependentAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "analysis_settings": {
                "module"  : "",
                "type"    : "",
                "settings": {}
            },
            "remove_modelparts": true,
            "remove_analysis"  : true
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.__remove_modelparts = parameters["remove_modelparts"].GetBool()
        self.__remove_analysis = parameters["remove_analysis"].GetBool()
        self.__existing_root_model_part_names = []
        self.__analysis = None

    def Initialize(self, _: dict):
        pass

    def InitializeIteration(self, _: dict):
        pass

    def Execute(self, _: dict):
        # make a record of existing root model parts
        self.__existing_root_model_part_names = [model_part_name for model_part_name in self.model.GetModelPartNames()]

        self.__analysis = RetrieveObject(self.model, self.parameters["analysis_settings"].Clone())
        self.__analysis.Run()

    def FinalizeIteration(self, _: dict):
        # delete the current analysis
        if self.__remove_analysis and self.__analysis is not None:
            del self.__analysis
            self.__analysis = None

        # remove all the modelparts created by the current execution
        if self.__remove_modelparts:
            for model_part_name in self.model.GetModelPartNames():
                if model_part_name not in self.__existing_root_model_part_names:
                    self.model.DeleteModelPart(model_part_name)

    def Finalize(self, _: dict):
        pass

    def GetAnalysis(self, _: dict):
        if self.__analysis is not None:
            return self.__analysis
        else:
            raise RuntimeError("Requesting analysis when it is not run for the current iteration.")

