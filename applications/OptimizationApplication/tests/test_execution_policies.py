
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator

class TestExecutionPolicies(kratos_unittest.TestCase):
    def test_IndependentAnalysisExecutionPolicy(self):
        model = Kratos.Model()
        parameters = Kratos.Parameters("""{
            "name"    : "test",
            "type"    : "independent_analysis_execution_policy",
            "settings": {
                "analysis_type"    : "orchestrators.SequentialOrchestrator",
                "analysis_settings": {
                    "orchestrator": {
                        "settings": {
                            "stage_checkpoints": false
                        }
                    },
                    "stages":[]
                }
            }
        }""")
        execution_policy = ExecutionPolicyDecorator(model, parameters, OptimizationProblem())
        execution_policy.Initialize()
        execution_policy.Execute()

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()