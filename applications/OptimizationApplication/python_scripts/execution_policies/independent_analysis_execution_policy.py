from typing import Union

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.multistage_analysis import MultistageAnalysis
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationProcessFactory

class IndependentAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__()

        self.model = model
        self.parameters = parameters

        default_settings = Kratos.Parameters("""{
            "analysis_module"  : "",
            "analysis_type"    : "",
            "analysis_settings": {}
        }""")
        self.parameters.ValidateAndAssignDefaults(default_settings)

    def Execute(self):
        current_analysis: Union[AnalysisStage, MultistageAnalysis] = OptimizationProcessFactory(self.parameters["analysis_module"].GetString(), self.parameters["analysis_type"].GetString(), self.model, self.parameters["analysis_settings"].Clone(), required_object_type=Union[AnalysisStage, MultistageAnalysis])
        current_analysis.Run()