import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import json
from KratosMultiphysics.OptimizationApplication.base_schema import ValidateJSON, BaseSchema

class TestOptimizationSchema(KratosUnittest.TestCase):
    def test_0(self):
        json_schema = BaseSchema.model_json_schema()
        with open('../python_scripts/schema/optimization_application_schema.json', 'w') as file:
            json.dump(json_schema, file, indent=4)  # Convert dictionary to JSON string and write to file

    def test_1(self):
        with open("algorithm_tests/analysis_based_tests/algorithm_gradient_projection/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

    def test_2(self):
        with open("algorithm_tests/analysis_based_tests/algorithm_nesterov_accelerated_gradient/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

    def test_3(self):
        with open("algorithm_tests/analysis_based_tests/algorithm_relaxed_gradient_projection/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

    def test_3(self):
        with open("algorithm_tests/analysis_based_tests/algorithm_relaxed_gradient_projection/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

    def test_4(self):
        with open("algorithm_tests/analysis_based_tests/algorithm_steepest_descent/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

    def test_5(self):
        with open("algorithm_tests/analysis_based_tests/algorithm_steepest_descent_qnbb/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

    def test_6(self):
        with open("algorithm_tests/nlopt_tests/mma_shell_thickness_opt/optimization_parameters.json", "r") as file_input:
            json_data = json.load(file_input)
            valid, validate_data = ValidateJSON(json_data)
            self.assertTrue(valid)

if __name__ == '__main__':
    KratosUnittest.main()