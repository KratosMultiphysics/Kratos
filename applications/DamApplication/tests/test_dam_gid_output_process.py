import json
import os
import re
import sys
import tempfile
from unittest import mock

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication  # noqa: F401
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis
from KratosMultiphysics.DamApplication.gid_dam_output_process import GiDDamOutputProcess
from KratosMultiphysics.gid_output_process import GiDOutputProcess

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)

MECHANICAL_CASE_DIR = os.path.join(TESTS_DIR, "joint_elastic_cohesive_2d_normal")
MECHANICAL_CASE = os.path.join(MECHANICAL_CASE_DIR, "joint_elastic_cohesive_2d_normal")

_RESULT_PATTERN = re.compile(r'^Result "([^"]+)" "[^"]+" (\S+) (\S+) (\S+)')


class TestDamGiDOutputProcess(KratosUnittest.TestCase):
    """Tests for the legacy GiD output adapter over the standard GiD writer."""

    def _load_parameters(self, problem_name, output_configuration,
                         start_time=0.0, end_time=2.0, time_step=1.0,
                         time_scale="Seconds"):
        with open(MECHANICAL_CASE + "_parameters.json", "r") as parameter_file:
            parameters = json.load(parameter_file)

        parameters["problem_data"]["problem_name"] = problem_name
        parameters["problem_data"]["start_time"] = start_time
        parameters["problem_data"]["end_time"] = end_time
        parameters["problem_data"]["time_step"] = time_step
        parameters["problem_data"]["time_scale"] = time_scale
        parameters["solver_settings"]["model_import_settings"]["input_filename"] = MECHANICAL_CASE
        parameters["output_configuration"] = output_configuration
        parameters["output_device_list"] = []
        parameters.pop("output_processes", None)
        return parameters

    def _run(self, parameters, working_directory):
        previous_directory = os.getcwd()
        os.chdir(working_directory)
        try:
            simulation = DamAnalysis(KM.Model(), KM.Parameters(json.dumps(parameters)))
            simulation.Run()
            return simulation
        finally:
            os.chdir(previous_directory)

    def _legacy_configuration(self, output_control_type="time_s", output_interval=1.0,
                              start_output_results=0.0, multi_file_flag="SingleFile",
                              gauss_point_results=None, skin_output=False,
                              plane_output=None, post_mode="GiD_PostAscii",
                              nodal_results=None, nodal_nonhistorical_results=None):
        return {
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": post_mode,
                    "WriteDeformedMeshFlag": "WriteDeformed",
                    "WriteConditionsFlag": "WriteElementsOnly",
                    "MultiFileFlag": multi_file_flag
                },
                "output_control_type": output_control_type,
                "output_interval": output_interval,
                "start_output_results": start_output_results,
                "body_output": True,
                "node_output": False,
                "skin_output": skin_output,
                "plane_output": plane_output or [],
                "nodal_results": (
                    ["DISPLACEMENT"] if nodal_results is None else nodal_results),
                "nodal_nonhistorical_results": (
                    [] if nodal_nonhistorical_results is None
                    else nodal_nonhistorical_results),
                "gauss_point_results": gauss_point_results or []
            },
            "point_data_configuration": []
        }

    def _result_entries(self, result_file):
        entries = []
        with open(result_file, "r") as file:
            for line in file:
                match = _RESULT_PATTERN.match(line.strip())
                if match:
                    entries.append(match.groups())
        return entries

    def _result_tags(self, result_file, variable_name="DISPLACEMENT"):
        return [
            tag for name, tag, _, _ in self._result_entries(result_file)
            if name == variable_name
        ]

    def _values_after_header(self, result_file, header):
        with open(result_file, "r") as file:
            lines = [line.strip() for line in file]

        header_index = lines.index(header)
        self.assertEqual(lines[header_index + 1], "Values")

        values = {}
        for line in lines[header_index + 2:]:
            if line == "End Values":
                break
            parts = line.split()
            if parts:
                values[int(parts[0])] = [float(value) for value in parts[1:]]
        return values

    def _list_files(self, working_directory):
        return sorted(
            name for name in os.listdir(working_directory)
            if name.endswith(".post.lst")
        )

    def _make_model_part(self):
        model = KM.Model()
        model_part = model.CreateModelPart("main")
        model_part.ProcessInfo.SetValue(KM.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KM.STEP, 0)
        return model_part

    def test_translation_of_legacy_time_units(self):
        cases = [
            ("time_s", 1.0, 86400.0),
            ("time_h", 3600.0, 24.0),
            ("time_d", 86400.0, 1.0),
            ("time_w", 604800.0, 1.0 / 7.0)
        ]

        for output_control_type, expected_factor, expected_start in cases:
            with self.subTest(output_control_type=output_control_type):
                with tempfile.TemporaryDirectory() as working_directory:
                    model_part = self._make_model_part()
                    import KratosMultiphysics.DamApplication as KratosDam
                    model_part.ProcessInfo.SetValue(
                        KratosDam.TIME_UNIT_CONVERTER, 86400.0)

                    process = GiDDamOutputProcess(
                        model_part,
                        os.path.join(working_directory, "translated"),
                        1.0,
                        KM.Parameters(json.dumps(self._legacy_configuration(
                            output_control_type=output_control_type,
                            output_interval=2.0,
                            start_output_results=5.0
                        )))
                    )

                    result_configuration = process.param["result_file_configuration"]
                    self.assertEqual(
                        result_configuration["output_control_type"].GetString(), "time")
                    self.assertAlmostEqual(
                        result_configuration["output_interval"].GetDouble(), 2.0 * expected_factor)
                    self.assertFalse(result_configuration.Has("start_output_results"))
                    self.assertEqual(
                        result_configuration["file_label"].GetString(), "time")
                    self.assertEqual(
                        result_configuration["time_label_format"].GetString(), "{:.12f}")
                    self.assertAlmostEqual(process._dam_output_factor, expected_factor)
                    self.assertAlmostEqual(process._dam_output_interval, 2.0)
                    self.assertAlmostEqual(process._dam_start_output_results, 5.0)
                    self.assertAlmostEqual(process._dam_start_time_output, expected_start)

    def test_empty_or_missing_configuration_uses_legacy_defaults(self):
        cases = [
            ("empty", KM.Parameters("{}")),
            ("missing", None)
        ]
        for name, settings in cases:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as working_directory:
                    process = GiDDamOutputProcess(
                        self._make_model_part(),
                        os.path.join(working_directory, name),
                        0.0,
                        settings
                    )

                    result_configuration = process.param["result_file_configuration"]
                    self.assertEqual(
                        result_configuration["output_control_type"].GetString(), "time")
                    self.assertAlmostEqual(
                        result_configuration["output_interval"].GetDouble(), 1.0)
                    self.assertAlmostEqual(process._dam_output_factor, 1.0)
                    self.assertAlmostEqual(process._dam_start_output_results, 0.0)
                    self.assertAlmostEqual(process._dam_start_time_output, 0.0)

    def test_output_frequency_alias_is_translated(self):
        with tempfile.TemporaryDirectory() as working_directory:
            configuration = self._legacy_configuration(output_interval=1.0)
            result_configuration = configuration["result_file_configuration"]
            del result_configuration["output_interval"]
            result_configuration["output_frequency"] = 3.0

            process = GiDDamOutputProcess(
                self._make_model_part(),
                os.path.join(working_directory, "frequency"),
                0.0,
                KM.Parameters(json.dumps(configuration))
            )

            self.assertFalse(
                process.param["result_file_configuration"].Has("output_frequency"))
            self.assertAlmostEqual(
                process.param["result_file_configuration"]["output_interval"].GetDouble(),
                3.0)
            self.assertAlmostEqual(process._dam_output_interval, 3.0)

    def test_step_control_is_translated_as_time_in_seconds(self):
        with tempfile.TemporaryDirectory() as working_directory:
            with mock.patch.object(KM.Logger, "PrintWarning") as warning:
                process = GiDDamOutputProcess(
                    self._make_model_part(),
                    os.path.join(working_directory, "step"),
                    0.0,
                    KM.Parameters(json.dumps(self._legacy_configuration(
                        output_control_type="step", output_interval=2.0)))
                )

            self.assertEqual(
                process.param["result_file_configuration"]["output_control_type"].GetString(),
                "time")
            self.assertAlmostEqual(process._dam_output_factor, 1.0)
            self.assertAlmostEqual(process._dam_output_interval, 2.0)
            warning.assert_called()

    def test_unknown_control_type_is_rejected(self):
        with tempfile.TemporaryDirectory() as working_directory:
            with self.assertRaises(Exception):
                GiDDamOutputProcess(
                    self._make_model_part(),
                    os.path.join(working_directory, "invalid"),
                    0.0,
                    KM.Parameters(json.dumps(self._legacy_configuration(
                        output_control_type="time_years")))
                )

    def test_output_time_labels_for_legacy_units(self):
        cases = [
            ("time_s", 1.0, 3.0, 1.0, ["0", "1", "2", "3"]),
            ("time_h", 1.0, 10800.0, 3600.0, ["0", "1", "2", "3"]),
            ("time_d", 1.0, 259200.0, 86400.0, ["0", "1", "2", "3"]),
            ("time_w", 1.0, 1209600.0, 604800.0, ["0", "1", "2"])
        ]

        for output_control_type, interval, end_time, time_step, expected_tags in cases:
            with self.subTest(output_control_type=output_control_type):
                with tempfile.TemporaryDirectory() as working_directory:
                    parameters = self._load_parameters(
                        "labels_" + output_control_type,
                        self._legacy_configuration(
                            output_control_type=output_control_type,
                            output_interval=interval),
                        end_time=end_time,
                        time_step=time_step)
                    simulation = self._run(parameters, working_directory)

                    self.assertIsInstance(simulation.gid_output, GiDOutputProcess)
                    self.assertEqual(
                        self._result_tags(os.path.join(
                            working_directory,
                            "labels_{0}.post.res".format(output_control_type))),
                        expected_tags)

    def test_fractional_output_labels_preserve_precision(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters(
                "fractional_labels",
                self._legacy_configuration(output_interval=1.0 / 3.0),
                end_time=1.0,
                time_step=1.0 / 3.0)
            self._run(parameters, working_directory)

            self.assertEqual(
                self._result_tags(os.path.join(
                    working_directory, "fractional_labels.post.res")),
                ["0", "0.3333333333333333", "0.6666666666666666", "1"])

    def test_scheduling_when_interval_and_time_step_differ(self):
        cases = [
            ("coarse_step", 2.0, 1.0, 6.0, ["0", "2", "4", "6"]),
            ("fine_step", 0.5, 2.0, 5.0, ["0", "2", "4"]),
            ("skipped_intervals", 5.0, 2.0, 10.0, ["0", "5", "10"])
        ]

        for name, time_step, interval, end_time, expected_tags in cases:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as working_directory:
                    parameters = self._load_parameters(
                        name,
                        self._legacy_configuration(output_interval=interval),
                        end_time=end_time,
                        time_step=time_step)
                    self._run(parameters, working_directory)

                    self.assertEqual(
                        self._result_tags(os.path.join(
                            working_directory, name + ".post.res")),
                        expected_tags)

    def test_delayed_and_nonzero_start_scheduling(self):
        cases = [
            ("delayed", 0.0, 3.0, 1.0, 1.0, ["1", "2", "3"]),
            ("nonzero_start", 1.5, 4.0, 0.5, 0.0, ["1.5", "2.5", "3.5"])
        ]

        for name, start_time, end_time, time_step, start_output, expected_tags in cases:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as working_directory:
                    parameters = self._load_parameters(
                        name,
                        self._legacy_configuration(
                            start_output_results=start_output),
                        start_time=start_time,
                        end_time=end_time,
                        time_step=time_step)
                    self._run(parameters, working_directory)

                    self.assertEqual(
                        self._result_tags(os.path.join(
                            working_directory, name + ".post.res")),
                        expected_tags)

    def test_multiple_files_initial_and_delayed_output(self):
        cases = [
            ("multi_initial", 0.0, ["0", "1", "2"]),
            ("multi_delayed", 1.0, ["1", "2"])
        ]

        for name, start_output, expected_tags in cases:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as working_directory:
                    parameters = self._load_parameters(
                        name,
                        self._legacy_configuration(
                            start_output_results=start_output,
                            multi_file_flag="MultipleFiles"),
                        end_time=2.0,
                        time_step=1.0)
                    self._run(parameters, working_directory)

                    for tag in expected_tags:
                        self.assertTrue(os.path.exists(os.path.join(
                            working_directory,
                            "{0}_{1}.post.res".format(name, tag))))
                    if "0" not in expected_tags:
                        self.assertFalse(os.path.exists(os.path.join(
                            working_directory, name + "_0.post.res")))

                    list_files = self._list_files(working_directory)
                    self.assertEqual(len(list_files), 1)
                    with open(os.path.join(working_directory, list_files[0]), "r") as file:
                        list_content = file.read()
                    self.assertIn("Multiple\n", list_content)
                    for tag in expected_tags:
                        self.assertIn(
                            "{0}_{1}.post.res".format(name, tag), list_content)

    def test_binary_single_and_multiple_files(self):
        for multi_file_flag, expected_files in [
            ("SingleFile", ["binary_single.post.bin"]),
            ("MultipleFiles", [
                "binary_multi_0.post.bin",
                "binary_multi_1.post.bin"
            ])
        ]:
            with self.subTest(multi_file_flag=multi_file_flag):
                with tempfile.TemporaryDirectory() as working_directory:
                    problem_name = (
                        "binary_single" if multi_file_flag == "SingleFile"
                        else "binary_multi"
                    )
                    parameters = self._load_parameters(
                        problem_name,
                        self._legacy_configuration(
                            multi_file_flag=multi_file_flag,
                            post_mode="GiD_PostBinary"),
                        end_time=1.0,
                        time_step=1.0)
                    self._run(parameters, working_directory)

                    for file_name in expected_files:
                        self.assertTrue(os.path.exists(os.path.join(
                            working_directory, file_name)))
                    self.assertEqual(len(self._list_files(working_directory)), 1)

    def test_step_control_multifile_uses_legacy_time_semantics(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters(
                "step_multi",
                self._legacy_configuration(
                    output_control_type="step",
                    multi_file_flag="MultipleFiles"),
                end_time=2.0,
                time_step=1.0)
            self._run(parameters, working_directory)

            for tag in ("0", "1", "2"):
                self.assertTrue(os.path.exists(os.path.join(
                    working_directory, "step_multi_{0}.post.res".format(tag))))
            self.assertEqual(len(self._list_files(working_directory)), 1)

    def test_mixed_simulation_and_output_units(self):
        cases = [
            ("mixed_hours", "time_h", 12.0, ["24", "36", "48", "60", "72", "84"]),
            ("matching_days", "time_d", 1.0, ["1", "2", "3"])
        ]

        for name, output_control_type, interval, expected_tags in cases:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as working_directory:
                    parameters = self._load_parameters(
                        name,
                        self._legacy_configuration(
                            output_control_type=output_control_type,
                            output_interval=interval),
                        start_time=1.0,
                        end_time=3.5,
                        time_step=0.5,
                        time_scale="Days")
                    self._run(parameters, working_directory)

                    self.assertEqual(
                        self._result_tags(os.path.join(
                            working_directory, name + ".post.res")),
                        expected_tags)

    def test_nodal_displacement_values(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters(
                "nodal_values",
                self._legacy_configuration(),
                end_time=2.0,
                time_step=1.0)
            self._run(parameters, working_directory)

            values = self._values_after_header(
                os.path.join(working_directory, "nodal_values.post.res"),
                'Result "DISPLACEMENT" "Kratos" 1 Vector OnNodes')
            self.assertAlmostEqual(values[4][0], 1.24864e-05, delta=1.0e-8)
            self.assertAlmostEqual(values[4][1], 8.74018e-06, delta=1.0e-8)
            self.assertAlmostEqual(values[4][2], 0.0, delta=1.0e-12)
            self.assertAlmostEqual(values[1][0], 0.0, delta=1.0e-12)

    def test_legacy_nodal_young_modulus_historical_request_is_translated(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters(
                "legacy_nodal_young_modulus",
                self._legacy_configuration(
                    nodal_results=["DISPLACEMENT", "NODAL_YOUNG_MODULUS"]),
                end_time=2.0,
                time_step=1.0)

            with mock.patch.object(KM.Logger, "PrintWarning") as warning:
                simulation = self._run(parameters, working_directory)

            result_file = os.path.join(
                working_directory, "legacy_nodal_young_modulus.post.res")
            self.assertEqual(
                self._result_tags(result_file, "YOUNG_MODULUS"), ["0", "1", "2"])
            self.assertEqual(
                self._result_tags(result_file, "NODAL_YOUNG_MODULUS"), [])
            self.assertEqual(
                self._result_tags(result_file, "DISPLACEMENT"), ["0", "1", "2"])

            for label in ("0", "1", "2"):
                values = self._values_after_header(
                    result_file,
                    'Result "YOUNG_MODULUS" "Kratos" {0} Scalar OnNodes'.format(label))
                self.assertEqual(len(values), 24)
                for value in values.values():
                    self.assertAlmostEqual(value[0], 0.0)

            result_configuration = simulation.gid_output.param["result_file_configuration"]
            nodal_results = result_configuration["nodal_results"].GetStringArray()
            self.assertIn("YOUNG_MODULUS", nodal_results)
            self.assertNotIn("NODAL_YOUNG_MODULUS", nodal_results)
            self.assertNotIn(
                "NODAL_YOUNG_MODULUS",
                result_configuration["nodal_nonhistorical_results"].GetStringArray())

            warning_messages = [str(call) for call in warning.call_args_list]
            self.assertTrue(any(
                "NODAL_YOUNG_MODULUS" in message for message in warning_messages))

    def test_legacy_nodal_young_modulus_nonhistorical_request_is_translated(self):
        cases = [
            ("nonhistorical_only", ["DISPLACEMENT"],
             ["NODAL_YOUNG_MODULUS"]),
            ("nonhistorical_with_existing_historical",
             ["DISPLACEMENT", "YOUNG_MODULUS"],
             ["NODAL_YOUNG_MODULUS"])
        ]

        for name, nodal_results, nonhistorical_results in cases:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as working_directory:
                    parameters = self._load_parameters(
                        name,
                        self._legacy_configuration(
                            nodal_results=nodal_results,
                            nodal_nonhistorical_results=nonhistorical_results),
                        end_time=2.0,
                        time_step=1.0)

                    with mock.patch.object(KM.Logger, "PrintWarning"):
                        simulation = self._run(parameters, working_directory)

                    result_file = os.path.join(working_directory, name + ".post.res")
                    self.assertEqual(
                        self._result_tags(result_file, "YOUNG_MODULUS"),
                        ["0", "1", "2"])
                    self.assertEqual(
                        self._result_tags(result_file, "NODAL_YOUNG_MODULUS"), [])

                    values = self._values_after_header(
                        result_file,
                        'Result "YOUNG_MODULUS" "Kratos" 0 Scalar OnNodes')
                    self.assertEqual(len(values), 24)
                    for value in values.values():
                        self.assertAlmostEqual(value[0], 0.0)

                    result_configuration = simulation.gid_output.param["result_file_configuration"]
                    translated_nodal_results = (
                        result_configuration["nodal_results"].GetStringArray())
                    self.assertEqual(
                        translated_nodal_results.count("YOUNG_MODULUS"), 1)
                    self.assertNotIn(
                        "NODAL_YOUNG_MODULUS", translated_nodal_results)
                    self.assertNotIn(
                        "NODAL_YOUNG_MODULUS",
                        result_configuration["nodal_nonhistorical_results"].GetStringArray())

    def test_nodal_young_modulus_translation_does_not_modify_input_parameters(self):
        with tempfile.TemporaryDirectory() as working_directory:
            input_parameters = KM.Parameters(json.dumps(self._legacy_configuration(
                nodal_results=["DISPLACEMENT", "NODAL_YOUNG_MODULUS"])))
            process = GiDDamOutputProcess(
                self._make_model_part(),
                os.path.join(working_directory, "input_parameters"),
                0.0,
                input_parameters)

            self.assertEqual(
                input_parameters["result_file_configuration"][
                    "nodal_results"].GetStringArray(),
                ["DISPLACEMENT", "NODAL_YOUNG_MODULUS"])
            self.assertEqual(
                process.param["result_file_configuration"][
                    "nodal_results"].GetStringArray(),
                ["DISPLACEMENT", "YOUNG_MODULUS"])

    def test_nodal_young_modulus_translation_handles_missing_nodal_results(self):
        configuration = self._legacy_configuration(
            nodal_nonhistorical_results=["NODAL_YOUNG_MODULUS"])
        del configuration["result_file_configuration"]["nodal_results"]

        with tempfile.TemporaryDirectory() as working_directory:
            process = GiDDamOutputProcess(
                self._make_model_part(),
                os.path.join(working_directory, "missing_nodal_results"),
                0.0,
                KM.Parameters(json.dumps(configuration)))

            result_configuration = process.param["result_file_configuration"]
            self.assertEqual(
                result_configuration["nodal_results"].GetStringArray(),
                ["YOUNG_MODULUS"])
            self.assertEqual(
                result_configuration["nodal_nonhistorical_results"].GetStringArray(),
                [])

    def test_gauss_point_output_is_delegated(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters(
                "gauss",
                self._legacy_configuration(
                    gauss_point_results=["CAUCHY_STRESS_VECTOR"]),
                end_time=2.0,
                time_step=1.0)
            self._run(parameters, working_directory)

            result_file = os.path.join(working_directory, "gauss.post.res")
            entries = self._result_entries(result_file)
            gauss_entries = [
                entry for entry in entries
                if entry[0] == "CAUCHY_STRESS_VECTOR" and entry[3] == "OnGaussPoints"
            ]
            self.assertEqual(sorted(set(entry[1] for entry in gauss_entries)), ["0", "1", "2"])
            with open(result_file, "r") as file:
                result_content = file.read()
            self.assertIn('GaussPoints "tri1_element_gp"', result_content)
            self.assertIn('GaussPoints "quad4_element_gp"', result_content)

            values = self._values_after_header(
                result_file,
                'Result "CAUCHY_STRESS_VECTOR" "Kratos" 1 Matrix OnGaussPoints "tri1_element_gp"')
            self.assertAlmostEqual(values[1][0], 0.0, delta=1.0e-6)
            self.assertAlmostEqual(values[1][1], -7726.33, delta=1.0)
            self.assertAlmostEqual(values[1][2], 5966.27, delta=1.0)

    def test_plane_output_is_delegated(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters(
                "plane",
                self._legacy_configuration(
                    plane_output=[{
                        "normal": [1.0, 0.0, 0.0],
                        "point": [1.0, 1.0, 0.0]
                    }]),
                end_time=2.0,
                time_step=1.0)
            self._run(parameters, working_directory)

            cut_result_file = os.path.join(working_directory, "plane_cuts.post.res")
            self.assertTrue(os.path.exists(cut_result_file))
            self.assertEqual(self._result_tags(cut_result_file), ["0", "1", "2"])


if __name__ == "__main__":
    KratosUnittest.main()
