
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper

class TestExecutionPolicies(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()

    def test_IndependentAnalysisExecutionPolicy(self):
        parameters = Kratos.Parameters("""{
            "name"                     : "test",
            "execution_policy_settings": {
                "type"    : "IndependentAnalysisExecutionPolicy",
                "settings": {
                    "analysis_settings": {
                        "module"  : "KratosMultiphysics",
                        "type"    : "MultistageAnalysis",
                        "settings": {
                            "stages": [],
                            "execution_list":[]
                        }
                    }
                }
            }
        }""")
        execution_policy = ExecutionPolicyWrapper(self.model, parameters)
        execution_policy.Initialize({})
        execution_policy.Execute({})

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()