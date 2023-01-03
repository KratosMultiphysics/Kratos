
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.utilities.execution_policies.independent_analysis_execution_policy import IndependentAnalysisExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.execution_policies.stepping_analysis_execution_policy import SteppingAnalysisExecutionPolicy

class TestExecutionPolicies(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

    def test_IndependentAnalysisExecutionPolicy(self):
        parameters = Kratos.Parameters("""{
            "log_in_file": false,
            "settings"       : {
                "analysis_name": "KratosMultiphysics.multistage_analysis.MultistageAnalysis",
                "analysis_parameters": {
                    "stages": []
                }
            }
        }""")
        execution_policy = IndependentAnalysisExecutionPolicy(self.model, parameters)
        execution_policy.Initialize({})
        execution_policy.Run({})

    # @kratos_unittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    # def test_SteppingAnalysisExecutionPolicy(self):
    #     parameters = Kratos.Parameters("""{
    #         "log_in_file": false,
    #         "settings"       : {
    #             "analysis_name": "KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis.StructuralMechanicsAnalysis",
    #             "analysis_parameters": {
    #                 "problem_data": {
    #                     "echo_level": 0,
    #                     "parallel_type": "OpenMP"
    #                 },
    #                 "solver_settings": {
    #                     "model_part_name": "Structure",
    #                     "solver_type": "Dynamic",
    #                     "domain_size": 3,
    #                     "time_stepping":{

    #                     }
    #                 }
    #             }
    #         }
    #     }""")
    #     execution_policy = SteppingAnalysisExecutionPolicy(self.model, parameters)
    #     execution_policy.Initialize({})
    #     execution_policy.Run({})

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    kratos_unittest.main()