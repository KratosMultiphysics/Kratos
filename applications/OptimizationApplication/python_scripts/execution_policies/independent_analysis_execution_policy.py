#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  License:		 BSD License
#					 license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#

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
            }
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

    def Initialize(self, _: dict):
        pass

    def Execute(self, _: dict):
        current_analysis = RetrieveObject(self.model, self.parameters["analysis_settings"].Clone())
        current_analysis.Run()