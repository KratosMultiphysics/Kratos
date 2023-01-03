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
from KratosMultiphysics.OptimizationApplication.utilities.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.execution_policies.execution_policy import RetrieveClass

class IndependentAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "analysis_name"      : "PLEASE_PROVIDE_FULL_MODULE_PATH.CLASS_NAME_WITH_RUN_METHOD",
            "analysis_parameters": {}
        }""")
        self.parameters["settings"].ValidateAndAssignDefaults(default_settings)
        self.analysis_class = RetrieveClass(self.parameters["settings"]["analysis_name"].GetString())
        self.analysis_parameters = self.parameters["settings"]["analysis_parameters"]

    def Initialize(self, _: dict):
        pass

    def RunAnalysis(self, _: dict):
        current_analysis_parameters = self.analysis_parameters.Clone()
        current_analysis = self.analysis_class(self.model, current_analysis_parameters)
        current_analysis.Run()