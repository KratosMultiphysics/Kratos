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
from KratosMultiphysics.DamApplication.dam_young_modulus_output_process import (
    DamYoungModulusOutputProcess,
)
from KratosMultiphysics.DamApplication.gid_dam_output_process import GiDDamOutputProcess

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
MECHANICAL_CASE = os.path.join(
    TESTS_DIR, "joint_elastic_cohesive_2d_normal", "joint_elastic_cohesive_2d_normal")

_RESULT_PATTERN = re.compile(r'^Result "([^"]+)" "([^"]+)" (\S+) (\S+) (\S+)$')


def _CreateSolidModelPart():
    model = KM.Model()
    model_part = model.CreateModelPart("main")
    model_part.AddNodalSolutionStepVariable(KM.YOUNG_MODULUS)

    coordinates = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))
    for node_id, (x, y) in enumerate(coordinates, start=1):
        model_part.CreateNewNode(node_id, x, y, 0.0)

    properties_1 = model_part.CreateNewProperties(1)
    properties_1[KM.YOUNG_MODULUS] = 1.0e9
    properties_1[KM.POISSON_RATIO] = 0.2
    properties_1[KM.DENSITY] = 1000.0

    properties_2 = model_part.CreateNewProperties(2)
    properties_2[KM.YOUNG_MODULUS] = 2.0e9
    properties_2[KM.POISSON_RATIO] = 0.2
    properties_2[KM.DENSITY] = 1000.0

    model_part.CreateNewElement(
        "SmallDisplacementThermoMechanicElement2D3N", 1, [1, 2, 3], properties_1)
    model_part.CreateNewElement(
        "SmallDisplacementThermoMechanicElement2D3N", 2, [2, 4, 3], properties_2)

    return model, model_part


def _CreateProcess(model, model_part_name="main"):
    settings = KM.Parameters()
    settings.AddString("model_part_name", model_part_name)
    return DamYoungModulusOutputProcess(model, settings)


def _GidResultBlocks(result_file):
    blocks = []
    current_block = None
    reading_values = False

    with open(result_file, "r") as file:
        for line in file:
            line = line.strip()
            match = _RESULT_PATTERN.match(line)
            if match:
                current_block = {
                    "name": match.group(1),
                    "label": match.group(3),
                    "values": {}
                }
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


class TestDamYoungModulusOutputProcess(KratosUnittest.TestCase):
    """Tests for the nonhistorical NODAL_YOUNG_MODULUS preparation."""

    def test_constant_material_property(self):
        model, model_part = _CreateSolidModelPart()
        for properties in model_part.Properties:
            properties[KM.YOUNG_MODULUS] = 1.0e9

        _CreateProcess(model).ExecuteBeforeSolutionLoop()

        for node in model_part.Nodes:
            self.assertAlmostEqual(node[KratosDam.NODAL_YOUNG_MODULUS], 1.0e9)
            self.assertEqual(node.GetSolutionStepValue(KM.YOUNG_MODULUS), 0.0)

        for properties in model_part.Properties:
            self.assertEqual(properties[KM.YOUNG_MODULUS], 1.0e9)
        for element in model_part.Elements:
            self.assertEqual(element[KM.YOUNG_MODULUS], 0.0)

    def test_multiple_materials_shared_nodes_are_averaged(self):
        model, model_part = _CreateSolidModelPart()

        _CreateProcess(model).ExecuteBeforeSolutionLoop()

        # Nodes 2 and 3 are shared by both materials. The visualization
        # projection is the arithmetic mean, so they show 1.5e9 while the
        # material-interior nodes keep their own value.
        expected_values = {
            1: 1.0e9,
            2: 1.5e9,
            3: 1.5e9,
            4: 2.0e9,
        }
        for node_id, expected_value in expected_values.items():
            self.assertAlmostEqual(
                model_part.GetNode(node_id)[KratosDam.NODAL_YOUNG_MODULUS],
                expected_value)

    def test_joint_elements_do_not_contribute(self):
        model, model_part = _CreateSolidModelPart()
        for properties in model_part.Properties:
            properties[KM.YOUNG_MODULUS] = 1.0e9

        joint_properties = model_part.CreateNewProperties(3)
        joint_properties[KM.DENSITY] = 1000.0

        model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        model_part.CreateNewNode(7, 3.0, 0.0, 0.0)
        model_part.CreateNewNode(8, 3.0, 1.0, 0.0)
        # The joint shares node 3 with the solid elements and adds the
        # joint-only nodes 5 to 8.
        model_part.CreateNewElement(
            "SmallDisplacementInterfaceElement2D4N", 3, [3, 5, 6, 8], joint_properties)

        _CreateProcess(model).ExecuteBeforeSolutionLoop()

        # Shared nodes receive the solid value, not a value diluted by the
        # joint. Joints without YOUNG_MODULUS do not contribute.
        for node_id in range(1, 5):
            self.assertAlmostEqual(
                model_part.GetNode(node_id)[KratosDam.NODAL_YOUNG_MODULUS], 1.0e9)

        # Nodes without a material-property source are left untouched.
        for node_id in range(5, 9):
            self.assertFalse(
                model_part.GetNode(node_id).Has(KratosDam.NODAL_YOUNG_MODULUS))

        self.assertFalse(joint_properties.Has(KM.YOUNG_MODULUS))

    def test_accessor_backed_materials_are_skipped(self):
        model, model_part = _CreateSolidModelPart()
        model_part.GetNode(1).SetSolutionStepValue(KM.YOUNG_MODULUS, 5.0e8)
        KratosDam.AddYoungModulusDatabaseAccessor(model_part)

        _CreateProcess(model).ExecuteBeforeSolutionLoop()

        for node in model_part.Nodes:
            self.assertFalse(node.Has(KratosDam.NODAL_YOUNG_MODULUS))
        self.assertEqual(
            model_part.GetNode(1).GetSolutionStepValue(KM.YOUNG_MODULUS), 5.0e8)

    def test_repeated_preparation_is_idempotent(self):
        model, model_part = _CreateSolidModelPart()
        process = _CreateProcess(model)

        # DamAnalysis calls ExecuteBeforeSolutionLoop twice (before and after
        # the solver initialization) and ExecuteBeforeOutputStep every step.
        process.ExecuteBeforeSolutionLoop()
        process.ExecuteBeforeOutputStep()
        process.ExecuteBeforeOutputStep()

        self.assertAlmostEqual(
            model_part.GetNode(1)[KratosDam.NODAL_YOUNG_MODULUS], 1.0e9)
        self.assertAlmostEqual(
            model_part.GetNode(4)[KratosDam.NODAL_YOUNG_MODULUS], 2.0e9)


class TestDamYoungModulusOutputIntegration(KratosUnittest.TestCase):
    """End-to-end checks through DamAnalysis, GiD and VTK."""

    def _load_parameters(self, problem_name, end_time=1.0):
        with open(MECHANICAL_CASE + "_parameters.json", "r") as parameter_file:
            parameters = json.load(parameter_file)

        parameters["problem_data"]["problem_name"] = problem_name
        parameters["problem_data"]["end_time"] = end_time
        parameters["problem_data"]["time_step"] = 1.0
        parameters["solver_settings"]["model_import_settings"]["input_filename"] = MECHANICAL_CASE
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

    def _legacy_configuration(self, nodal_results):
        return {
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostAscii",
                    "WriteDeformedMeshFlag": "WriteDeformed",
                    "WriteConditionsFlag": "WriteElementsOnly",
                    "MultiFileFlag": "SingleFile"
                },
                "output_control_type": "time_s",
                "output_interval": 1.0,
                "start_output_results": 0,
                "body_output": True,
                "node_output": False,
                "skin_output": False,
                "plane_output": [],
                "nodal_results": nodal_results,
                "gauss_point_results": []
            },
            "point_data_configuration": []
        }

    def _vtk_process(self, output_path, nodal_data_value_variables):
        return {
            "python_module": "vtk_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "VtkOutputProcess",
            "Parameters": {
                "model_part_name": "MainModelPart.mechanical_computing_domain",
                "output_control_type": "time",
                "output_interval": 1.0,
                "file_format": "ascii",
                "output_precision": 7,
                "output_path": output_path,
                "save_output_files_in_folder": True,
                "write_deformed_configuration": True,
                "nodal_data_value_variables": nodal_data_value_variables
            }
        }

    def _gid_process(self, output_name):
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
                        "output_interval": 1.0,
                        "body_output": True,
                        "node_output": False,
                        "skin_output": False,
                        "plane_output": [],
                        "nodal_results": ["DISPLACEMENT"],
                        "nodal_nonhistorical_results": ["NODAL_YOUNG_MODULUS"],
                        "gauss_point_results": []
                    },
                    "point_data_configuration": []
                }
            }
        }

    def _vtk_values(self, vtk_file, variable_name):
        with open(vtk_file, "r") as file:
            lines = [line.strip() for line in file]

        header_index = next(
            index for index, line in enumerate(lines)
            if line.startswith(variable_name + " "))
        number_of_values = int(lines[header_index].split()[2])

        values = [
            float(line)
            for line in lines[header_index + 1: header_index + 1 + number_of_values]
        ]
        self.assertEqual(len(values), number_of_values)
        return values

    def test_legacy_gid_nodal_young_modulus(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("legacy_nodal_young_modulus")
            parameters["output_configuration"] = self._legacy_configuration(
                ["DISPLACEMENT", "NODAL_YOUNG_MODULUS"])

            with mock.patch.object(KM.Logger, "PrintWarning") as warning:
                self._run(parameters, working_directory)

            result_file = os.path.join(
                working_directory, "legacy_nodal_young_modulus.post.res")
            blocks = _GidResultBlocks(result_file)

            young_modulus_blocks = [
                block for block in blocks if block["name"] == "NODAL_YOUNG_MODULUS"]
            self.assertEqual(
                [float(block["label"]) for block in young_modulus_blocks],
                [0.0, 1.0])

            for block in young_modulus_blocks:
                self.assertEqual(sorted(block["values"]), list(range(1, 25)))
                for value in block["values"].values():
                    self.assertAlmostEqual(value, 1.0e9)

            self.assertTrue(any(
                block["name"] == "DISPLACEMENT" for block in blocks))

            warning_messages = [str(call) for call in warning.call_args_list]
            self.assertTrue(any(
                "NODAL_YOUNG_MODULUS" in message for message in warning_messages))

    def test_vtk_nodal_young_modulus(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("vtk_nodal_young_modulus", end_time=2.0)
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(
                    vtk_output_path, ["NODAL_YOUNG_MODULUS"])]
            }

            self._run(parameters, working_directory)

            vtk_files = sorted(
                os.path.join(vtk_output_path, file_name)
                for file_name in os.listdir(vtk_output_path)
                if file_name.endswith(".vtk")
            )
            self.assertEqual(len(vtk_files), 2)
            for vtk_file in vtk_files:
                values = self._vtk_values(vtk_file, "NODAL_YOUNG_MODULUS")
                self.assertEqual(len(values), 24)
                for value in values:
                    self.assertAlmostEqual(value, 1.0e9)

    def test_standard_gid_nonhistorical_nodal_young_modulus(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("standard_gid_nodal_young_modulus")
            parameters["output_processes"] = {
                "gid_output": [self._gid_process("standard_gid_nodal_young_modulus")]
            }

            simulation = self._run(parameters, working_directory)

            self.assertTrue(any(
                isinstance(process, DamYoungModulusOutputProcess)
                for process in simulation.list_of_processes))

            result_file = os.path.join(
                working_directory, "standard_gid_nodal_young_modulus.post.res")
            young_modulus_blocks = [
                block for block in _GidResultBlocks(result_file)
                if block["name"] == "NODAL_YOUNG_MODULUS"]
            self.assertTrue(young_modulus_blocks)
            for block in young_modulus_blocks:
                self.assertEqual(len(block["values"]), 24)
                for value in block["values"].values():
                    self.assertAlmostEqual(value, 1.0e9)

    def test_no_young_modulus_output_requested(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("no_young_modulus")
            parameters["output_configuration"] = self._legacy_configuration(
                ["DISPLACEMENT"])
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(vtk_output_path, ["DISPLACEMENT"])]
            }

            simulation = self._run(parameters, working_directory)

            self.assertFalse(any(
                isinstance(process, DamYoungModulusOutputProcess)
                for process in simulation.list_of_processes))

            result_file = os.path.join(
                working_directory, "no_young_modulus.post.res")
            self.assertFalse(any(
                block["name"] == "NODAL_YOUNG_MODULUS"
                for block in _GidResultBlocks(result_file)))

            for vtk_file in os.listdir(vtk_output_path):
                with open(os.path.join(vtk_output_path, vtk_file), "r") as file:
                    self.assertNotIn("NODAL_YOUNG_MODULUS", file.read())

    def test_explicit_process_is_not_duplicated(self):
        with tempfile.TemporaryDirectory() as working_directory:
            parameters = self._load_parameters("explicit_preparation")
            parameters["output_configuration"] = self._legacy_configuration(
                ["NODAL_YOUNG_MODULUS"])
            vtk_output_path = os.path.join(working_directory, "VTK_Output")
            parameters["output_processes"] = {
                "vtk_output": [self._vtk_process(
                    vtk_output_path, ["NODAL_YOUNG_MODULUS"])]
            }
            parameters["constraints_process_list"].append({
                "python_module": "dam_young_modulus_output_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "process_name": "DamYoungModulusOutputProcess",
                "Parameters": {
                    "model_part_name": "MainModelPart.mechanical_computing_domain"
                }
            })

            simulation = self._run(parameters, working_directory)

            preparation_processes = [
                process for process in simulation.list_of_processes
                if isinstance(process, DamYoungModulusOutputProcess)
            ]
            self.assertEqual(len(preparation_processes), 1)

            result_file = os.path.join(
                working_directory, "explicit_preparation.post.res")
            initial_block = next(
                block for block in _GidResultBlocks(result_file)
                if block["name"] == "NODAL_YOUNG_MODULUS" and block["label"] == "0")
            for value in initial_block["values"].values():
                self.assertAlmostEqual(value, 1.0e9)


if __name__ == "__main__":
    KratosUnittest.runTests(
        KratosUnittest.TestLoader().loadTestsFromTestCase(
            TestDamYoungModulusOutputProcess))
