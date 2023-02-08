
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator

class TestExecutionPolicies(kratos_unittest.TestCase):
    def test_IndependentAnalysisExecutionPolicy(self):
        model = Kratos.Model()
        parameters = Kratos.Parameters("""{
            "execution_policy_name"    : "test",
            "execution_policy_type"    : "IndependentAnalysisExecutionPolicy",
            "execution_policy_settings": {
                "analysis_module"  : "KratosMultiphysics",
                "analysis_type"    : "MultistageAnalysis",
                "analysis_settings": {
                    "stages": [],
                    "execution_list":[]
                }
            }
        }""")
        execution_policy = ExecutionPolicyDecorator(model, parameters)
        execution_policy.ExecuteInitialize()
        execution_policy.Execute()

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()