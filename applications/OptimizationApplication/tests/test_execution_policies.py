import os, shutil
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator

class TestExecutionPolicies(kratos_unittest.TestCase):
    def test_IndependentAnalysisExecutionPolicy(self):
        model = Kratos.Model()
        model.CreateModelPart("MainModelPart")

        parameters = Kratos.Parameters("""{
            "name"    : "test",
            "type"    : "independent_analysis_execution_policy",
            "settings": {
                "analysis_type"    : "orchestrators.SequentialOrchestrator",
                "analysis_model_part_name": "MainModelPart",
                "analysis_settings": {
                    "input_folder"                : "temp",
                    "output_folder"               : "auxiliary_files/<name>/<step>",
                    "project_parameters_file_name": "orchestrator_parameters.json",
                    "output_settings": {
                        "hdf5_file_name": "dummy.h5",
                        "list_of_values": [

                        ]
                    }
                }
            }
        }""")

        # Create a temporary folder and copy the orchestrator_parameters.json file into it
        if not os.path.exists("temp"):
            os.makedirs("temp")
        shutil.copyfile("auxiliary_files/orchestrator_parameters.json", "temp/orchestrator_parameters.json")

        execution_policy = ExecutionPolicyDecorator(model, parameters, OptimizationProblem())
        execution_policy.Initialize()
        execution_policy.Execute()

        # Delete the temporary folder
        shutil.rmtree("temp")
        shutil.rmtree("auxiliary_files/test")

if __name__ == "__main__":
    kratos_unittest.main()