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
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import RetrieveClass

class IndependentAnalysisExecutionPolicy:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "analysis_type"      : "PLEASE_PROVIDE_FULL_MODULE_PATH.CLASS_NAME_WITH_RUN_METHOD",
            "analysis_parameters": {}
        }""")
        self.model = model

        parameters.ValidateAndAssignDefaults(default_settings)
        self.analysis_class = RetrieveClass(parameters["analysis_type"].GetString())
        self.analysis_parameters = parameters["analysis_parameters"]

    def Initialize(self, _: dict):
        pass

    def Execute(self, _: dict):
        current_analysis_parameters = self.analysis_parameters.Clone()
        current_analysis = self.analysis_class(self.model, current_analysis_parameters)
        current_analysis.Run()