from pathlib import Path
import shutil
import tempfile

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem


class TestExecutionPolicies(kratos_unittest.TestCase):
    @kratos_unittest.skipIfApplicationsNotAvailable("HDF5Application")
    def test_IndependentAnalysisExecutionPolicy(self):
        model = Kratos.Model()
        model.CreateModelPart("MainModelPart")

        test_dir = Path(__file__).resolve().parent
        source_parameters = test_dir / "auxiliary_files" / "orchestrator_parameters.json"

        with tempfile.TemporaryDirectory(prefix="kratos_execution_policy_") as temp_dir:
            temp_dir = Path(temp_dir)
            input_dir = temp_dir / "input"
            output_dir = temp_dir / "output"
            input_dir.mkdir()

            shutil.copyfile(source_parameters, input_dir / "orchestrator_parameters.json")

            parameters = Kratos.Parameters(
                f"""
                {{
                    "name": "test",
                    "type": "independent_analysis_execution_policy",
                    "settings": {{
                        "analysis_type": "orchestrators.SequentialOrchestrator",
                        "analysis_model_part_name": "MainModelPart",
                        "analysis_settings": {{
                            "input_folder": "{input_dir.as_posix()}",
                            "output_folder": "{(output_dir / "<name>" / "<step>").as_posix()}",
                            "project_parameters_file_name": "orchestrator_parameters.json",
                            "output_settings": {{
                                "hdf5_file_name": "dummy.h5",
                                "list_of_values": []
                            }}
                        }}
                    }}
                }}
                """
            )

            execution_policy = ExecutionPolicyDecorator(model, parameters, OptimizationProblem())
            execution_policy.Initialize()
            execution_policy.Execute()

if __name__ == "__main__":
    kratos_unittest.main()