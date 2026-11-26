import json
import os
import sys
import tempfile

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)

import dam_recording_lifecycle_process as recording_lifecycle_process


BASE_CASE = os.path.join(
    TESTS_DIR,
    "joint_elastic_cohesive_2d_normal",
    "joint_elastic_cohesive_2d_normal",
)


class TestDamAnalysisProcessLifecycle(KratosUnittest.TestCase):

    def _load_parameters(self, problem_name):
        with open(BASE_CASE + "_parameters.json", "r") as parameter_file:
            parameters = json.load(parameter_file)

        parameters["problem_data"]["problem_name"] = problem_name
        parameters["problem_data"]["end_time"] = 2.0
        parameters["problem_data"]["time_step"] = 1.0
        parameters["problem_data"]["number_of_threads"] = 1
        parameters["solver_settings"]["echo_level"] = 0
        parameters["solver_settings"]["model_import_settings"][
            "input_filename"] = BASE_CASE

        parameters.pop("output_configuration", None)
        parameters.pop("output_processes", None)
        parameters["output_device_list"] = []
        return parameters

    def _run(self, parameters):
        with tempfile.TemporaryDirectory() as working_directory:
            previous_directory = os.getcwd()
            os.chdir(working_directory)
            try:
                model = KM.Model()
                simulation = DamAnalysis(model, KM.Parameters(json.dumps(parameters)))
                simulation.Run()
                return simulation
            finally:
                os.chdir(previous_directory)

    def test_ordinary_process_lifecycle_callbacks(self):
        parameters = self._load_parameters("lifecycle_callbacks")
        parameters["constraints_process_list"].insert(0, {
            "python_module": "dam_recording_lifecycle_process",
            "Parameters": {},
        })

        self._run(parameters)

        calls = recording_lifecycle_process.LAST_PROCESS.calls

        self.assertEqual(calls.count("ExecuteInitialize"), 1)
        self.assertEqual(calls.count("ExecuteBeforeSolutionLoop"), 1)
        self.assertEqual(calls.count("ExecuteInitializeSolutionStep"), 2)
        self.assertEqual(calls.count("ExecuteFinalizeSolutionStep"), 2)
        self.assertEqual(calls.count("ExecuteFinalize"), 1)

        self.assertEqual(calls[0], "ExecuteInitialize")
        self.assertLess(
            calls.index("ExecuteInitialize"),
            calls.index("ExecuteBeforeSolutionLoop"))
        self.assertLess(
            calls.index("ExecuteBeforeSolutionLoop"),
            calls.index("ExecuteInitializeSolutionStep"))
        self.assertLess(
            calls.index("ExecuteInitializeSolutionStep"),
            calls.index("ExecuteFinalizeSolutionStep"))
        self.assertEqual(calls[-1], "ExecuteFinalize")

    def test_standard_process_execute_initialize_is_applied(self):
        parameters = self._load_parameters("standard_process_initialize")
        parameters["constraints_process_list"].insert(0, {
            "python_module": "process_factory",
            "kratos_module": "KratosMultiphysics",
            "process_name": "ApplyConstantScalarValueProcess",
            "Parameters": {
                "model_part_name": "MainModelPart",
                "variable_name": "YOUNG_MODULUS",
                "value": 123.0,
                "is_fixed": False,
            },
        })

        simulation = self._run(parameters)

        computing_model_part = simulation.solver.GetComputingModelPart()
        self.assertGreater(computing_model_part.NumberOfNodes(), 0)
        for node in computing_model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(KM.YOUNG_MODULUS),
                123.0)


if __name__ == "__main__":
    KratosUnittest.main()
