import json
import os
import re
import sys
import tempfile
from unittest import mock

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis
from KratosMultiphysics.DamApplication.gid_dam_output_process import GiDDamOutputProcess
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)

MECHANICAL_CASE_DIR = os.path.join(TESTS_DIR, "joint_elastic_cohesive_2d_normal")
MECHANICAL_CASE = os.path.join(MECHANICAL_CASE_DIR, "joint_elastic_cohesive_2d_normal")
CONSTRUCTION_CASE_DIR = os.path.join(TESTS_DIR, "construction")
CONSTRUCTION_CASE = os.path.join(CONSTRUCTION_CASE_DIR, "construction")


class TestDamOutputProcesses(KratosUnittest.TestCase):
    """Tests for the standard output_processes integration in DamAnalysis."""

    def _load_parameters(self, problem_name, base_name=MECHANICAL_CASE):
        with open(base_name + "_parameters.json", "r") as parameter_file:
            parameters = json.load(parameter_file)

        parameters["problem_data"]["problem_name"] = problem_name
        parameters["problem_data"]["end_time"] = 2.0
        parameters["problem_data"]["time_step"] = 1.0
        parameters["solver_settings"]["model_import_settings"]["input_filename"] = base_name

        parameters.pop("output_configuration", None)
        parameters.pop("output_processes", None)
        parameters["output_device_list"] = []
        return parameters

    def _run(self, parameters, working_directory):
        previous_directory = os.getcwd()
        os.chdir(working_directory)
        try:
            model = KM.Model()
            simulation = DamAnalysis(model, KM.Parameters(json.dumps(parameters)))
            simulation.Run()
            return simulation
        finally:
            os.chdir(previous_directory)

    def _legacy_configuration(self, multi_file="SingleFile"):
        return {
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostAscii",
                    "WriteDeformedMeshFlag": "WriteDeformed",
                    "WriteConditionsFlag": "WriteElementsOnly",
                    "MultiFileFlag": multi_file
                },
                "output_control_type": "time_s",
                "output_interval": 1.0,
                "start_output_results": 0,
                "body_output": True,
                "node_output": False,
                "skin_output": False,
                "plane_output": [],
                "nodal_results": ["DISPLACEMENT"],
                "gauss_point_results": []
            },
            "point_data_configuration": []
        }

    def _gid_process(self, output_name, interval=1.0):
        return {
            "python_module": "gid_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "GiDOutputProcess",
            "Parameters": {
                "model_part_name": "MainModelPart.mechanical_computing_domain",
                "output_name": output_name,
                "postprocess_parameters": {
                    "result_file_configuration": {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostAscii",
                            "WriteDeformedMeshFlag": "WriteDeformed",
                            "WriteConditionsFlag": "WriteElementsOnly",
                            "MultiFileFlag": "SingleFile"
                        },
                        "file_label": "time",
                        "time_label_format": "{:.12f}",
                        "output_control_type": "time",
                        "output_interval": interval,
                        "body_output": True,
                        "node_output": False,
                        "skin_output": False,
                        "plane_output": [],
                        "nodal_results": ["DISPLACEMENT"],
                        "gauss_point_results": []
                    },
                    "point_data_configuration": []
                }
            }
        }

    def _vtk_process(self, output_path, variables=None, interval=1.0):
        if variables is None:
            variables = ["DISPLACEMENT"]
        return {
            "python_module": "vtk_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "VtkOutputProcess",
            "Parameters": {
                "model_part_name": "MainModelPart.mechanical_computing_domain",
                "output_control_type": "time",
                "output_interval": interval,
                "file_format": "ascii",
                "output_precision": 7,
                "output_path": output_path,
                "save_output_files_in_folder": True,
                "write_deformed_configuration": True,
                "nodal_solution_step_data_variables": variables
            }
        }

    def _result_tags(self, result_file, variable_name="DISPLACEMENT"):
        tags = []
        pattern = re.compile(
            r'^Result "' + re.escape(variable_name) +
            r'" "Kratos" (\S+) (?:Scalar|Vector|Matrix) OnNodes$')
        with open(result_file, "r") as file:
            for line in file:
                match = pattern.match(line.strip())
                if match:
                    tags.append(match.group(1))
        return tags

    def _scalar_result_blocks(self, result_file, variable_name):
        pattern = re.compile(
            r'^Result "' + re.escape(variable_name) +
            r'" "[^"]+" (\S+) Scalar OnNodes$')
        blocks = []
        current_block = None
        reading_values = False

        with open(result_file, "r") as file:
            for line in file:
                line = line.strip()
                match = pattern.match(line)
                if match:
                    current_block = {"label": match.group(1), "values": {}}
                    blocks.append(current_block)
                    reading_values = False
                    continue
                if current_block is None:
                    continue
                if line == "Values":
                    reading_values = True
                    continue
                if line == "End Values":
                    current_block = None
                    reading_values = False
                    continue
                if reading_values and line:
                    parts = line.split()
                    current_block["values"][int(parts[0])] = float(parts[1])
        return blocks

    def _vtk_values(self, vtk_file, variable_name):
        with open(vtk_file, "r") as file:
            lines = [line.strip() for line in file]

        header_index = next(
            index for index, line in enumerate(lines)
            if line.startswith(variable_name + " "))
        number_of_values = int(lines[header_index].split()[2])
        return [
            float(line)
            for line in lines[header_index + 1: header_index + 1 + number_of_values]
        ]

    def _vtk_files(self, working_directory):
        vtk_directory = os.path.join(working_directory, "VTK_Output")
        if not os.path.isdir(vtk_directory):
            return []
        return sorted(
            os.path.join(vtk_directory, file_name)
            for file_name in os.listdir(vtk_directory)
            if file_name.endswith(".vtk")
        )

    def _assert_vtk_field(self, vtk_file, variable_name):
        with open(vtk_file, "r") as file:
            lines = file.read().splitlines()

        self.assertTrue(lines[0].startswith("# vtk DataFile Version"))
        self.assertIn("ASCII", lines)
        self.assertTrue(any(line.startswith("POINT_DATA ") for line in lines))

        field_header = None
        for line in lines:
            if line.startswith(variable_name + " "):
                field_header = line
                break
        self.assertIsNotNone(field_header, "Missing VTK field " + variable_name)
        header_parts = field_header.split()
        number_of_components = int(header_parts[1])
        number_of_values = int(header_parts[2])

        field_index = lines.index(field_header)
        values = []
        for line in lines[field_index + 1:]:
            parts = line.split()
            if len(parts) != number_of_components:
                break
            values.extend(float(value) for value in parts)
        self.assertEqual(len(values), number_of_components * number_of_values)
        self.assertTrue(any(abs(value) > 0.0 for value in values))

    def _assert_recording_lifecycle(self, calls, output_steps):
        expected_counts = {
            "ExecuteInitialize": 1,
            "ExecuteBeforeSolutionLoop": 1,
            "ExecuteInitializeSolutionStep": output_steps,
            "ExecuteFinalizeSolutionStep": output_steps,
            "ExecuteBeforeOutputStep": output_steps,
            "PrintOutput": output_steps,
            "ExecuteAfterOutputStep": output_steps,
            "ExecuteFinalize": 1
        }
        for call_name, expected_count in expected_counts.items():
            self.assertEqual(calls.count(call_name), expected_count)
        self.assertEqual(calls[0], "ExecuteInitialize")
        self.assertEqual(calls[-1], "ExecuteFinalize")

        for index, call_name in enumerate(calls):
            if call_name == "PrintOutput":
                self.assertEqual(calls[index - 1], "ExecuteBeforeOutputStep")
                self.assertEqual(calls[index + 1], "ExecuteAfterOutputStep")

    def test_legacy_gid_only(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_gid")
            parameters["output_configuration"] = self._legacy_configuration()

            simulation = self._run(parameters, working_directory)

            self.assertIsNotNone(simulation.gid_output)
            self.assertEqual(simulation.list_of_standard_output_processes, [])
            self.assertEqual(
                self._result_tags(os.path.join(working_directory, "legacy_gid.post.res")),
                ["0", "1", "2"]
            )

    def test_empty_legacy_output_configuration(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("empty_legacy")
            parameters["output_configuration"] = {}

            simulation = self._run(parameters, working_directory)

            self.assertIsNotNone(simulation.gid_output)
            self.assertEqual(simulation.list_of_standard_output_processes, [])
            self.assertTrue(os.path.exists(
                os.path.join(working_directory, "empty_legacy.post.bin")))

    def test_legacy_gid_multiple_files(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_multi")
            parameters["output_configuration"] = self._legacy_configuration("MultipleFiles")

            simulation = self._run(parameters, working_directory)

            self.assertIsNotNone(simulation.gid_output)
            for label in ("0", "1", "2"):
                self.assertTrue(os.path.exists(
                    os.path.join(working_directory, "legacy_multi_{0}.post.res".format(label))))
            list_files = [
                name for name in os.listdir(working_directory)
                if name.endswith(".post.lst")
            ]
            self.assertEqual(len(list_files), 1)
            with open(os.path.join(working_directory, list_files[0]), "r") as list_file:
                list_content = list_file.read()
            self.assertIn("Multiple\n", list_content)
            for label in ("0", "1", "2"):
                self.assertIn(
                    "legacy_multi_{0}.post.res".format(label), list_content)

    def test_legacy_gid_cleans_previous_post_files(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_cleanup")
            parameters["output_configuration"] = self._legacy_configuration()

            stale_files = [
                "stale.post.res",
                "stale.post.msh",
                "stale.post.lst",
                "info.time"
            ]
            for file_name in stale_files:
                with open(os.path.join(working_directory, file_name), "w") as file:
                    file.write("stale")

            simulation = self._run(parameters, working_directory)

            self.assertIsNotNone(simulation.gid_output)
            for file_name in stale_files:
                self.assertFalse(os.path.exists(os.path.join(working_directory, file_name)))
            self.assertTrue(os.path.exists(
                os.path.join(working_directory, "legacy_cleanup.post.res")))

    def test_vtk_only(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("vtk_only")
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(vtk_output_path)]
            }

            sentinel_file = os.path.join(working_directory, "existing.post.res")
            with open(sentinel_file, "w") as file:
                file.write("existing GiD results")

            simulation = self._run(parameters, working_directory)

            self.assertIsNone(simulation.gid_output)
            self.assertEqual(len(simulation.list_of_standard_output_processes), 1)
            self.assertIsInstance(simulation.list_of_standard_output_processes[0], VtkOutputProcess)
            self.assertTrue(os.path.exists(sentinel_file))
            self.assertFalse(os.path.exists(os.path.join(working_directory, "vtk_only.post.res")))

            vtk_files = self._vtk_files(working_directory)
            self.assertEqual(
                [os.path.basename(vtk_file) for vtk_file in vtk_files],
                [
                    "mechanical_computing_domain_0_1.0000000.vtk",
                    "mechanical_computing_domain_0_2.0000000.vtk"
                ]
            )
            self._assert_vtk_field(vtk_files[0], "DISPLACEMENT")
            self._assert_vtk_field(vtk_files[1], "DISPLACEMENT")

    def test_legacy_gid_and_standard_vtk(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_and_vtk")
            parameters["output_configuration"] = self._legacy_configuration()
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "gid_output": [],
                "vtk_output": [self._vtk_process(vtk_output_path)]
            }

            simulation = self._run(parameters, working_directory)

            self.assertIsNotNone(simulation.gid_output)
            self.assertEqual(len(simulation.list_of_standard_output_processes), 1)
            self.assertIsInstance(simulation.list_of_standard_output_processes[0], VtkOutputProcess)
            self.assertEqual(
                self._result_tags(os.path.join(working_directory, "legacy_and_vtk.post.res")),
                ["0", "1", "2"]
            )
            self.assertEqual(len(self._vtk_files(working_directory)), 2)

    def test_independent_gid_and_vtk_scheduling(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("independent_schedule")
            parameters["problem_data"]["start_time"] = 1.5
            parameters["problem_data"]["time_step"] = 0.5
            parameters["problem_data"]["end_time"] = 4.0
            parameters["output_processes"] = {
                "gid_output": [self._gid_process("scheduled_gid", interval=2.0)],
                "vtk_output": [self._vtk_process(
                    os.path.join(working_directory, "VTK_Output"),
                    interval=3.0
                )]
            }

            simulation = self._run(parameters, working_directory)

            self.assertIsNone(simulation.gid_output)
            self.assertEqual(
                self._result_tags(os.path.join(working_directory, "scheduled_gid.post.res")),
                ["2", "4"]
            )
            self.assertEqual(
                [os.path.basename(vtk_file) for vtk_file in self._vtk_files(working_directory)],
                ["mechanical_computing_domain_0_3.0000000.vtk"]
            )
            self.assertFalse(os.path.exists(
                os.path.join(working_directory, "independent_schedule.post.res")))

    def test_standard_gid_and_standard_vtk(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_name")
            parameters["output_processes"] = {
                "gid_output": [self._gid_process("standard_gid")],
                "vtk_output": [self._vtk_process(os.path.join(working_directory, "VTK_Output"))]
            }

            simulation = self._run(parameters, working_directory)

            self.assertIsNone(simulation.gid_output)
            self.assertEqual(len(simulation.list_of_standard_output_processes), 2)
            self.assertIsInstance(simulation.list_of_standard_output_processes[0], GiDOutputProcess)
            self.assertIsInstance(simulation.list_of_standard_output_processes[1], VtkOutputProcess)
            self.assertEqual(
                self._result_tags(os.path.join(working_directory, "standard_gid.post.res")),
                ["1", "2"]
            )
            self.assertFalse(os.path.exists(os.path.join(working_directory, "legacy_name.post.res")))
            self.assertEqual(len(self._vtk_files(working_directory)), 2)

    def test_standard_gid_takes_precedence_over_legacy(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_name")
            parameters["output_configuration"] = self._legacy_configuration()
            parameters["output_processes"] = {
                "gid_output": [self._gid_process("standard_gid")]
            }

            with mock.patch.object(KM.Logger, "PrintWarning") as warning:
                simulation = self._run(parameters, working_directory)

            self.assertIsNone(simulation.gid_output)
            self.assertTrue(os.path.exists(os.path.join(working_directory, "standard_gid.post.res")))
            self.assertFalse(os.path.exists(os.path.join(working_directory, "legacy_name.post.res")))
            warning.assert_called()
            warning_messages = [str(call) for call in warning.call_args_list]
            self.assertTrue(any("Both 'output_configuration'" in message for message in warning_messages))

    def test_non_gid_process_named_gid_output_keeps_legacy_gid(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("fake_gid_group")
            parameters["output_configuration"] = self._legacy_configuration()
            parameters["output_processes"] = {
                "gid_output": [{
                    "python_module": "dam_recording_output_process",
                    "Parameters": {"record_name": "not_gid"}
                }]
            }

            with mock.patch.object(KM.Logger, "PrintWarning") as warning:
                simulation = self._run(parameters, working_directory)

            self.assertIsInstance(simulation.gid_output, GiDDamOutputProcess)
            self.assertEqual(len(simulation.list_of_standard_output_processes), 1)
            self.assertEqual(
                self._result_tags(os.path.join(working_directory, "fake_gid_group.post.res")),
                ["0", "1", "2"]
            )
            warning_messages = [str(call) for call in warning.call_args_list]
            self.assertFalse(any(
                "Both 'output_configuration'" in message for message in warning_messages))

    def test_missing_output_blocks(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("no_output")
            parameters.pop("output_device_list", None)

            simulation = self._run(parameters, working_directory)

            self.assertIsNone(simulation.gid_output)
            self.assertEqual(simulation.list_of_standard_output_processes, [])
            self.assertEqual(
                [name for name in os.listdir(working_directory) if ".post." in name],
                []
            )

    def test_monitoring_devices(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("devices")
            parameters["output_device_list"] = [{
                "python_module": "point_output_process",
                "kratos_module": "KratosMultiphysics",
                "process_name": "PointOutputProcess",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "position": [0.75, 0.25, 0.0],
                    "output_variables": ["DISPLACEMENT"],
                    "output_file_settings": {
                        "file_name": "monitor",
                        "output_path": "MonitoringDevices"
                    }
                }
            }]

            simulation = self._run(parameters, working_directory)

            self.assertIsNone(simulation.gid_output)
            monitor_file = os.path.join(working_directory, "MonitoringDevices", "monitor.dat")
            self.assertTrue(os.path.exists(monitor_file))
            with open(monitor_file, "r") as file:
                monitor_content = file.read()
            self.assertIn("DISPLACEMENT_X", monitor_content)
            self.assertIn("\n1.0 ", monitor_content)
            self.assertIn("\n2.0 ", monitor_content)

    def test_standard_output_process_lifecycle(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("lifecycle")
            parameters["output_processes"] = {
                "recording": [{
                    "python_module": "dam_recording_output_process",
                    "Parameters": {"record_name": "lifecycle"}
                }]
            }

            simulation = self._run(parameters, working_directory)

            self.assertEqual(len(simulation.list_of_standard_output_processes), 1)
            self._assert_recording_lifecycle(
                simulation.list_of_standard_output_processes[0].calls, 2)

    def test_multiple_output_processes_in_same_group(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("multiple_group")
            parameters["output_processes"] = {
                "recording": [
                    {
                        "python_module": "dam_recording_output_process",
                        "Parameters": {"record_name": "first"}
                    },
                    {
                        "python_module": "dam_recording_output_process",
                        "Parameters": {"record_name": "second"}
                    }
                ]
            }

            simulation = self._run(parameters, working_directory)

            self.assertEqual(len(simulation.list_of_standard_output_processes), 2)
            self.assertEqual(
                [process.settings["record_name"].GetString()
                 for process in simulation.list_of_standard_output_processes],
                ["first", "second"]
            )
            for process in simulation.list_of_standard_output_processes:
                self._assert_recording_lifecycle(process.calls, 2)

    def test_vtk_temperature_output(self):
        with tempfile.TemporaryDirectory() as working_directory:
            with open(CONSTRUCTION_CASE + "_parameters.json", "r") as parameter_file:
                parameters = json.load(parameter_file)

            parameters["problem_data"]["problem_name"] = "temperature_vtk"
            parameters["problem_data"]["end_time"] = 7200.0
            parameters["problem_data"]["time_step"] = 3600.0
            parameters["solver_settings"]["model_import_settings"]["input_filename"] = CONSTRUCTION_CASE
            parameters["construction_process"]["construction_input_file_name"] = os.path.join(
                CONSTRUCTION_CASE_DIR, "construction.txt")
            parameters["construction_process"]["ambient_input_file_name"] = os.path.join(
                CONSTRUCTION_CASE_DIR, "ambient.txt")
            parameters.pop("output_configuration", None)
            parameters.pop("output_device_list", None)
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(
                    vtk_output_path,
                    variables=["DISPLACEMENT", "TEMPERATURE"],
                    interval=3600.0
                )]
            }

            self._run(parameters, working_directory)

            vtk_files = self._vtk_files(working_directory)
            self.assertEqual(
                [os.path.basename(vtk_file) for vtk_file in vtk_files],
                [
                    "mechanical_computing_domain_0_3600.0000000.vtk",
                    "mechanical_computing_domain_0_7200.0000000.vtk"
                ]
            )
            self._assert_vtk_field(vtk_files[0], "DISPLACEMENT")
            self._assert_vtk_field(vtk_files[0], "TEMPERATURE")
            self._assert_vtk_field(vtk_files[1], "DISPLACEMENT")
            self._assert_vtk_field(vtk_files[1], "TEMPERATURE")

    def test_property_based_young_modulus_is_not_projected(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("property_young_modulus")
            parameters["output_configuration"] = self._legacy_configuration()
            parameters["output_configuration"]["result_file_configuration"][
                "nodal_results"] = ["DISPLACEMENT", "YOUNG_MODULUS"]
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(
                    vtk_output_path, variables=["YOUNG_MODULUS"])]
            }

            simulation = self._run(parameters, working_directory)

            self.assertIsInstance(simulation.gid_output, GiDDamOutputProcess)
            self.assertFalse(any(
                "YoungModulusOutputProcess" in type(process).__name__
                for process in simulation.list_of_processes))

            computing_model_part = simulation.solver.GetComputingModelPart()
            for node in computing_model_part.Nodes:
                self.assertEqual(node.GetSolutionStepValue(KM.YOUNG_MODULUS), 0.0)
                self.assertFalse(node.Has(KratosDam.NODAL_YOUNG_MODULUS))
            for properties in computing_model_part.Properties:
                if properties.Has(KM.YOUNG_MODULUS):
                    self.assertAlmostEqual(properties[KM.YOUNG_MODULUS], 1.0e9)

            result_file = os.path.join(
                working_directory, "property_young_modulus.post.res")
            self.assertEqual(
                self._result_tags(result_file, "YOUNG_MODULUS"), ["0", "1", "2"])
            self.assertEqual(
                self._result_tags(result_file, "DISPLACEMENT"), ["0", "1", "2"])
            for block in self._scalar_result_blocks(result_file, "YOUNG_MODULUS"):
                self.assertEqual(len(block["values"]), 24)
                for value in block["values"].values():
                    self.assertEqual(value, 0.0)

            vtk_files = self._vtk_files(working_directory)
            self.assertEqual(len(vtk_files), 2)
            for vtk_file in vtk_files:
                values = self._vtk_values(vtk_file, "YOUNG_MODULUS")
                self.assertEqual(len(values), 24)
                for value in values:
                    self.assertEqual(value, 0.0)

    def test_genuine_nodal_young_modulus_output(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("nodal_young_modulus")
            parameters["problem_data"]["end_time"] = 1.0
            parameters["solver_settings"]["mechanical_solver_settings"][
                "move_mesh_flag"] = False
            parameters["constraints_process_list"].insert(0, {
                "python_module": "impose_nodal_young_modulus_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "process_name": "ImposeNodalYoungModulusProcess",
                "Parameters": {
                    "model_part_name": "MainModelPart.mechanical_computing_domain",
                    "variable_name": "NODAL_YOUNG_MODULUS",
                    "Young_Modulus_1": 1.0e9,
                    "Young_Modulus_2": 1.0e9,
                    "Young_Modulus_3": 0.0,
                    "Young_Modulus_4": 0.0
                }
            })
            parameters["output_configuration"] = self._legacy_configuration()
            parameters["output_configuration"]["result_file_configuration"][
                "nodal_results"] = ["DISPLACEMENT", "YOUNG_MODULUS"]
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(
                    vtk_output_path, variables=["YOUNG_MODULUS"])]
            }

            simulation = self._run(parameters, working_directory)

            self.assertIsInstance(simulation.gid_output, GiDDamOutputProcess)
            self.assertFalse(any(
                "YoungModulusOutputProcess" in type(process).__name__
                for process in simulation.list_of_processes))

            computing_model_part = simulation.solver.GetComputingModelPart()
            expected_values = {
                node.Id: 1.0e9 + 1.0e9 * node.X
                for node in computing_model_part.Nodes
            }
            for node in computing_model_part.Nodes:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.YOUNG_MODULUS),
                    expected_values[node.Id])
                self.assertFalse(node.Has(KratosDam.NODAL_YOUNG_MODULUS))
            for properties in computing_model_part.Properties:
                if properties.Has(KM.YOUNG_MODULUS):
                    self.assertTrue(properties.HasAccessor(KM.YOUNG_MODULUS))
                    self.assertAlmostEqual(properties[KM.YOUNG_MODULUS], 1.0e9)

            result_file = os.path.join(
                working_directory, "nodal_young_modulus.post.res")
            young_modulus_blocks = self._scalar_result_blocks(
                result_file, "YOUNG_MODULUS")
            self.assertEqual(
                [float(block["label"]) for block in young_modulus_blocks],
                [0.0, 1.0])
            for block in young_modulus_blocks:
                self.assertEqual(
                    sorted(block["values"]), sorted(expected_values))
                for node_id, expected_value in expected_values.items():
                    self.assertAlmostEqual(
                        block["values"][node_id], expected_value)

            vtk_files = self._vtk_files(working_directory)
            self.assertEqual(len(vtk_files), 1)
            vtk_values = self._vtk_values(vtk_files[0], "YOUNG_MODULUS")
            self.assertEqual(len(vtk_values), len(expected_values))
            for value, node in zip(vtk_values, computing_model_part.Nodes):
                self.assertAlmostEqual(value, expected_values[node.Id])
