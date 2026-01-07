import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics.json_output_process import JsonOutputProcess

import itertools
import json
import os
from pathlib import Path

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestJsonOutputProcess(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("Main")

    def tearDown(self):
        KratosUtils.DeleteFileIfExisting(GetFilePath("external_forces_output.json"))
        KratosUtils.DeleteFileIfExisting(GetFilePath("resultant_external_forces_output.json"))

    def add_three_nodes_to_model_part(self):
        node_ids_and_coordinates = [
            (1, 0.0, 0.0, 0.0),
            (2, 1.0, 0.0, 0.0),
            (3, 1.0, 1.0, 0.0),
        ]
        for node_id, x, y, z in node_ids_and_coordinates:
            self.model_part.CreateNewNode(node_id, x, y, z)

    def add_vectors_to_nodes_of_model_part(self, vector_variable, nodal_vectors):
        for node, vector in zip(self.model_part.Nodes, nodal_vectors):
            node.SetSolutionStepValue(vector_variable, vector)

    def produce_json_output_file(self, process_settings):
        output_process = JsonOutputProcess(self.model, process_settings)

        output_process.ExecuteInitialize()

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 1.0)
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.5)

        output_process.ExecuteBeforeSolutionLoop()
        output_process.ExecuteFinalizeSolutionStep()

    def test_output_of_vector_variable(self):
        vector_variable_for_testing = KratosMultiphysics.EXTERNAL_FORCES_VECTOR
        self.model_part.AddNodalSolutionStepVariable(vector_variable_for_testing)
        self.add_three_nodes_to_model_part()
        test_vectors = [[1.0] * 3, [2.0] * 3, [3.0] * 3]
        self.add_vectors_to_nodes_of_model_part(
            vector_variable_for_testing, test_vectors
        )

        output_file_path = (
            Path(__file__).resolve().parent / "external_forces_output.json"
        )
        process_settings = KratosMultiphysics.Parameters(
            f"""{{
            "model_part_name": "{self.model_part.Name}",
            "output_file_name": "{output_file_path.as_posix()}",
            "historical_value": true,
            "output_variables": ["{vector_variable_for_testing.Name()}"]
        }}"""
        )
        self.produce_json_output_file(process_settings)

        with open(output_file_path, "r") as output_file:
            output_data = json.load(output_file)

        for node, vector in zip(self.model_part.Nodes, test_vectors):
            vector_variable_results_for_all_time_steps = output_data[f"NODE_{node.Id}"][
                vector_variable_for_testing.Name()
            ]
            self.assertEqual(len(vector_variable_results_for_all_time_steps), 1)
            self.assertVectorAlmostEqual(
                vector_variable_results_for_all_time_steps[0], vector
            )

    def test_output_of_resultant_vector_variable(self):
        vector_variable_for_testing = KratosMultiphysics.EXTERNAL_FORCES_VECTOR
        self.model_part.AddNodalSolutionStepVariable(vector_variable_for_testing)
        self.add_three_nodes_to_model_part()
        test_vectors = [[1.0] * 3, [2.0] * 3, [3.0] * 3]
        self.add_vectors_to_nodes_of_model_part(
            vector_variable_for_testing, test_vectors
        )
        ref_solution = [6.0, 6.0, 6.0]

        output_file_path = (
            Path(__file__).resolve().parent / "resultant_external_forces_output.json"
        )
        process_settings = KratosMultiphysics.Parameters(
            f"""{{
            "model_part_name": "{self.model_part.Name}",
            "output_file_name": "{output_file_path.as_posix()}",
            "historical_value": true,
            "output_variables": ["{vector_variable_for_testing.Name()}"],
            "resultant_solution": true
        }}"""
        )
        self.produce_json_output_file(process_settings)

        with open(output_file_path, "r") as output_file:
            output_data = json.load(output_file)

        resultant_vector_variable_results_for_all_time_steps = output_data["RESULTANT"][
            vector_variable_for_testing.Name()
        ]
        self.assertEqual(len(resultant_vector_variable_results_for_all_time_steps), 1)
        self.assertVectorAlmostEqual(
            resultant_vector_variable_results_for_all_time_steps[0],
            ref_solution,
        )


if __name__ == "__main__":
    KratosUnittest.main()
