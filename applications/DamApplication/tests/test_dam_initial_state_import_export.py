"""Lightweight selective initial-state import/export tests."""
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.restart_utility import RestartUtility

# Make the source-tree Dam scripts visible when running from the build tree.
Dam.__path__.insert(0, str(Path(__file__).resolve().parents[1] / "python_scripts"))

from KratosMultiphysics.DamApplication import dam_initial_state_file
from KratosMultiphysics.DamApplication.save_dam_initial_state_process import (
    SaveDamInitialStateProcess)
from KratosMultiphysics.DamApplication.load_dam_initial_state_process import (
    LoadDamInitialStateProcess)
import test_dam_selective_initialization as workflow_tests
import test_dam_selective_initial_stress_geometries as geometry_tests


class TestDamInitialStateImportExport(KratosUnittest.TestCase):

    @staticmethod
    def _save(model, file_name, stresses=True, temperature=False, reference=False,
              stress_source="stored_initial_stress", model_part="MainModelPart"):
        settings = KM.Parameters(json.dumps({
            "model_part_name": model_part,
            "save_initial_stresses": stresses,
            "save_initial_temperature": temperature,
            "save_reference_temperature": reference,
            "stress_source": stress_source,
            "output_file_name": str(file_name),
            "echo_level": 0}))
        SaveDamInitialStateProcess(model, settings).ExecuteFinalize()

    @staticmethod
    def _load(model, file_name, stresses=True, temperature=False, reference=False,
              policy="preserve", model_part="MainModelPart"):
        settings = KM.Parameters(json.dumps({
            "model_part_name": model_part,
            "input_file_name": str(file_name),
            "transfer_initial_stresses": stresses,
            "transfer_initial_temperature": temperature,
            "transfer_reference_temperature": reference,
            "reference_temperature_policy": policy,
            "echo_level": 0}))
        LoadDamInitialStateProcess(model, settings).ExecuteBeforeSolutionLoop()

    @staticmethod
    def _max_error(actual, expected):
        if isinstance(expected, list):
            return max((TestDamInitialStateImportExport._max_error(a, e)
                        for a, e in zip(actual, expected)), default=0.0)
        return abs(actual - expected)

    @staticmethod
    def _stress_map_error(actual, expected):
        return max(
            abs(a - b)
            for element_id in actual
            for actual_point, expected_point in zip(
                actual[element_id], expected[element_id])
            for a, b in zip(actual_point, expected_point))

    @staticmethod
    def _add_save(settings, file_name, stress_source, stresses=True,
                  temperature=False, reference=False):
        settings["loads_process_list"].append({
            "python_module": "save_dam_initial_state_process",
            "kratos_module": "KratosMultiphysics.DamApplication",
            "Parameters": {
                "model_part_name": "MainModelPart.Parts_Concrete",
                "save_initial_stresses": stresses,
                "save_initial_temperature": temperature,
                "save_reference_temperature": reference,
                "stress_source": stress_source,
                "output_file_name": str(file_name)}})

    @staticmethod
    def _add_load(settings, file_name, stresses, temperature, reference, policy):
        settings["loads_process_list"].append({
            "python_module": "load_dam_initial_state_process",
            "kratos_module": "KratosMultiphysics.DamApplication",
            "Parameters": {
                "model_part_name": "MainModelPart.Parts_Concrete",
                "input_file_name": str(file_name),
                "transfer_initial_stresses": stresses,
                "transfer_initial_temperature": temperature,
                "transfer_reference_temperature": reference,
                "reference_temperature_policy": policy}})

    def _run_full(self, directory, settings, thermal, probe_stored=True,
                  expected_error=None, temperature_offset=0.0):
        """Writes a case and executes a real, unmodified DamAnalysis.Run()."""
        directory.mkdir(parents=True, exist_ok=True)
        if thermal:
            mdpa = workflow_tests._MDPA.replace(
                "LinearElasticPlaneStrain2DLaw", "ThermalLinearElastic2DPlaneStrain")
            mdpa = mdpa.replace(
                "THICKNESS 1.0",
                "THICKNESS 1.0\nTHERMAL_EXPANSION 0.00001\n"
                "CONDUCTIVITY 1.0\nSPECIFIC_HEAT 1.0")
            mdpa += ("\nBegin Elements EulerianConvDiff2D4N\n"
                     "3 1 1 2 4 3\n4 1 3 4 6 5\nEnd Elements\n")
            mdpa += ("\nBegin SubModelPart Thermal_Concrete\nBegin SubModelPartNodes\n"
                     "1\n2\n3\n4\n5\n6\nEnd SubModelPartNodes\nBegin SubModelPartElements\n"
                     "3\n4\nEnd SubModelPartElements\nEnd SubModelPart\n")
            mdpa += ("\nBegin SubModelPart Thermal_Boundary\nBegin SubModelPartNodes\n"
                     "1\n2\n5\n6\nEnd SubModelPartNodes\nEnd SubModelPart\n")
            mdpa += "\nBegin NodalData TEMPERATURE\n"
            mdpa += "".join("{} 0 {}\n".format(
                i, 20 + 2 * ((i - 1) // 2) + temperature_offset)
                for i in range(1, 7))
            mdpa += "End NodalData\n"
            for variable, value in (("DENSITY", 2.), ("SPECIFIC_HEAT", 1.),
                                    ("CONDUCTIVITY", 1.)):
                mdpa += "\nBegin NodalData " + variable + "\n"
                mdpa += "".join("{} 0 {}\n".format(i, value) for i in range(1, 7))
                mdpa += "End NodalData\n"
            (directory / "reference_temperature.txt").write_text(
                "".join("{} {}\n".format(i, 20 + 2 * ((i - 1) // 2))
                        for i in range(1, 7)))
        else:
            mdpa = workflow_tests._MDPA
        (directory / "mini.mdpa").write_text(mdpa)

        probe_file = directory / "probe.json"
        settings["loads_process_list"].append({
            "python_module": "state_probe_process",
            "Parameters": {
                "model_part_name": "MainModelPart.Parts_Concrete",
                "output_file_name": str(probe_file),
                "probe_stored_initial_stress": probe_stored}})
        (directory / "parameters.json").write_text(json.dumps(settings))

        runner = (Path(__file__).resolve().parent / "selective_initial_stress" /
                  "run_full.py")
        source_scripts = Path(__file__).resolve().parents[1] / "python_scripts"
        result = subprocess.run(
            [sys.executable, str(runner), str(source_scripts), "parameters.json"],
            cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=300, env=os.environ.copy())
        (directory / "run.log").write_text(result.stdout)
        if expected_error is not None:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(expected_error, result.stdout)
            return None
        self.assertEqual(result.returncode, 0, result.stdout)
        return json.loads(probe_file.read_text())

    def test_file_roundtrip_2d_3d(self):
        helper = geometry_tests.TestDamSelectiveInitialStressGeometries()
        for name in ("Q4", "H8"):
            with self.subTest(geometry=name):
                geometry = geometry_tests.GEOMETRIES[name]
                dim = geometry["dim"]
                law = helper._law(dim, thermal=True)
                source_model, source, source_nodes, _ = helper._build(geometry, law)
                helper._set_quadratic_displacement(
                    source_nodes, geometry["coords"], dim)
                export_model, export_part, export_nodes, export_element = helper._build(
                    geometry, law)
                Dam.SelectiveInitialStressUtility.AssignInitialStresses(
                    source, export_part)
                temperature = {node.Id: 20.0 + node.Id for node in export_nodes}
                for node in export_nodes:
                    node.SetSolutionStepValue(KM.TEMPERATURE, temperature[node.Id])
                    node.SetSolutionStepValue(
                        Dam.NODAL_REFERENCE_TEMPERATURE, temperature[node.Id])
                exported_stress = helper._stresses(export_element, export_part)

                with tempfile.TemporaryDirectory(prefix="dam_state_roundtrip_") as tmp:
                    path = Path(tmp) / "state.json"
                    self._save(export_model, path, True, True, True)
                    data = json.loads(path.read_text())
                    self.assertEqual(
                        data["fields"],
                        ["initial_stresses", "temperature", "reference_temperature"])

                    import_model, import_part, import_nodes, import_element = helper._build(
                        geometry, law)
                    self._load(import_model, path, True, True, True)
                    error = self._max_error(
                        helper._stresses(import_element, import_part), exported_stress)
                    self.assertLess(error, 1.e-12)
                    for node in import_nodes:
                        self.assertAlmostEqual(
                            node.GetSolutionStepValue(KM.TEMPERATURE),
                            temperature[node.Id], delta=1.e-14)
                        self.assertAlmostEqual(
                            node.GetSolutionStepValue(
                                Dam.NODAL_REFERENCE_TEMPERATURE),
                            temperature[node.Id], delta=1.e-14)
                        for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                            self.assertEqual(value, 0.0)

    def test_independent_field_selection(self):
        helper = geometry_tests.TestDamSelectiveInitialStressGeometries()
        geometry = geometry_tests.GEOMETRIES["Q4"]
        law = helper._law(2)
        source_model, source, source_nodes, _ = helper._build(geometry, law)
        helper._set_quadratic_displacement(source_nodes, geometry["coords"], 2)
        export_model, export_part, export_nodes, export_element = helper._build(
            geometry, law)
        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, export_part)
        for node in export_nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, 20.0 + node.Id)
            node.SetSolutionStepValue(
                Dam.NODAL_REFERENCE_TEMPERATURE, 30.0 + node.Id)

        with tempfile.TemporaryDirectory(prefix="dam_state_fields_") as tmp:
            tmp = Path(tmp)
            stress_file = tmp / "stress.json"
            self._save(export_model, stress_file, True, False, False, "current_cauchy")
            import_model, import_part, import_nodes, import_element = helper._build(
                geometry, law)
            self._load(import_model, stress_file, True, False, False)
            self.assertLess(self._max_error(
                helper._stresses(import_element, import_part),
                helper._stresses(export_element, export_part)), 1.e-12)

            temperature_file = tmp / "temperature.json"
            self._save(export_model, temperature_file, False, True, False)
            import_model, import_part, import_nodes, _ = helper._build(geometry, law)
            self._load(import_model, temperature_file, False, True, False)
            for node in import_nodes:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KM.TEMPERATURE),
                    20.0 + node.Id, delta=1.e-14)

            reference_file = tmp / "reference.json"
            self._save(export_model, reference_file, False, False, True)
            import_model, import_part, import_nodes, _ = helper._build(geometry, law)
            self._load(import_model, reference_file, False, False, True)
            for node in import_nodes:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE),
                    30.0 + node.Id, delta=1.e-14)

    def test_stored_initial_stress_persistence(self):
        helper = geometry_tests.TestDamSelectiveInitialStressGeometries()
        geometry = geometry_tests.GEOMETRIES["Q4"]
        law = helper._law(2, thermal=True)
        source_model, source, source_nodes, _ = helper._build(geometry, law)
        helper._set_quadratic_displacement(source_nodes, geometry["coords"], 2)
        model, part, nodes, element = helper._build(geometry, law)
        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, part)
        stored_before = Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(part)
        current_before = helper._stresses(element, part)

        for node in nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, 30.0)
        stored_after = Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(part)
        current_after = helper._stresses(element, part)
        self.assertEqual(
            self._stress_map_error(stored_after, stored_before), 0.0)
        self.assertGreater(
            self._max_error(current_after, current_before), 1.e-6)

        with tempfile.TemporaryDirectory(prefix="dam_state_persistence_") as tmp:
            tmp = Path(tmp)
            path = tmp / "stored.json"
            self._save(model, path, True, False, False, "stored_initial_stress")
            current_path = tmp / "current.json"
            self._save(model, current_path, True, False, False, "current_cauchy")
            current_data = json.loads(current_path.read_text())
            self.assertLess(self._stress_map_error(
                current_data["initial_stresses"],
                {str(element.Id): [list(value) for value in values]
                 for element, values in
                 [(element, helper._stresses(element, part)) for element in part.Elements]}),
                1.e-12)
            self.assertGreater(
                self._stress_map_error(
                    current_data["initial_stresses"],
                    {str(element.Id): [list(value) for value in values]
                     for element, values in
                     [(element, stored_before[element.Id]) for element in part.Elements]}),
                1.e-6)
            import_model, import_part, _, import_element = helper._build(geometry, law)
            self._load(import_model, path, True, False, False)
            self.assertEqual(
                self._stress_map_error(
                    Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(
                        import_part), stored_before), 0.0)

            part.ProcessInfo[KM.TIME] = 1.0
            part.ProcessInfo[KM.STEP] = 1
            restart_settings = KM.Parameters(json.dumps({
                "input_filename": str(tmp / "MainModelPart"),
                "input_output_path": "checkpoints",
                "restart_save_frequency": 0.0,
                "save_restart_files_in_folder": True,
                "restart_control_type": "time",
                "restart_load_file_label": "1.0",
                "load_restart_files_from_folder": True,
                "echo_level": 0}))
            RestartUtility(part, restart_settings).SaveRestart()
            loaded_model, loaded_part, _, _ = helper._build(geometry, law)
            RestartUtility(loaded_part, restart_settings).LoadRestart()
            self.assertEqual(
                self._stress_map_error(
                    Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(
                        loaded_part), stored_before), 0.0)

    def test_reference_temperature_policies(self):
        helper = geometry_tests.TestDamSelectiveInitialStressGeometries()
        geometry = geometry_tests.GEOMETRIES["Q4"]
        law = helper._law(2, thermal=True)
        source_model, source, source_nodes, _ = helper._build(geometry, law)
        helper._set_quadratic_displacement(source_nodes, geometry["coords"], 2)
        export_model, export_part, export_nodes, _ = helper._build(geometry, law)
        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, export_part)
        for node in export_nodes:
            node.SetSolutionStepValue(KM.TEMPERATURE, 25.0)
            node.SetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE, 20.0)

        with tempfile.TemporaryDirectory(prefix="dam_state_policy_") as tmp:
            path = Path(tmp) / "state.json"
            self._save(export_model, path, True, True, True, "stored_initial_stress")

            import_model, import_part, import_nodes, import_element = helper._build(
                geometry, law)
            self._load(import_model, path, True, True, True, "preserve")
            for node in import_nodes:
                self.assertEqual(node.GetSolutionStepValue(KM.TEMPERATURE), 25.0)
                self.assertEqual(
                    node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE), 20.0)
                for value in node.GetSolutionStepValue(KM.DISPLACEMENT):
                    self.assertEqual(value, 0.0)
            imported = Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(
                import_part)
            actual = {
                str(import_element.Id):
                helper._stresses(import_element, import_part)}
            self.assertEqual(
                self._stress_map_error(
                    actual,
                    {str(key): value for key, value in imported.items()}),
                0.0)

            import_model, import_part, import_nodes, _ = helper._build(geometry, law)
            self._load(import_model, path, False, True, True, "preserve")
            for node in import_nodes:
                self.assertEqual(node.GetSolutionStepValue(KM.TEMPERATURE), 25.0)
                self.assertEqual(
                    node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE), 20.0)

            import_model, import_part, import_nodes, import_element = helper._build(
                geometry, law)
            self._load(import_model, path, True, True, False, "initial_temperature")
            for node in import_nodes:
                self.assertEqual(node.GetSolutionStepValue(KM.TEMPERATURE), 25.0)
                self.assertEqual(
                    node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE), 25.0)
            self.assertLess(self._stress_map_error(
                Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(import_part),
                Dam.SelectiveInitialStressUtility.GetStoredInitialStresses(export_part)),
                1.e-12)

            import_model, import_part, _, _ = helper._build(geometry, law)
            with self.assertRaisesRegex(RuntimeError, "Cannot both import"):
                self._load(import_model, path, True, True, True, "initial_temperature")
            with self.assertRaisesRegex(RuntimeError, "requires"):
                self._load(import_model, path, True, False, False, "initial_temperature")

    def test_negative_file_validation(self):
        helper = geometry_tests.TestDamSelectiveInitialStressGeometries()
        geometry = geometry_tests.GEOMETRIES["Q4"]
        law = helper._law(2)
        source_model, source, source_nodes, _ = helper._build(geometry, law)
        helper._set_quadratic_displacement(source_nodes, geometry["coords"], 2)
        export_model, export_part, _, _ = helper._build(geometry, law)
        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, export_part)

        with tempfile.TemporaryDirectory(prefix="dam_state_negative_") as tmp:
            tmp = Path(tmp)
            path = tmp / "state.json"
            self._save(export_model, path, True, False, False, "stored_initial_stress")

            def rejected(model, pattern, **kwargs):
                with self.assertRaisesRegex(RuntimeError, pattern):
                    self._load(model, path, **kwargs)

            model, part, nodes, element = helper._build(geometry, law)
            nodes[3].X0 += 0.1
            rejected(model, "reference coordinate")

            model, part, _, element = helper._build(
                geometry, law, connectivity=[2, 3, 4, 1])
            rejected(model, "connectivity")

            model, part, _, element = helper._build(
                geometry, law, integration_order=1)
            rejected(model, "integration rule")

            model, part, _, element = helper._build(
                geometry, helper._law(2, plane_stress=True))
            rejected(model, "constitutive law")

            model, part, _, element = helper._build(geometry, law)
            part.Properties[1][KM.YOUNG_MODULUS] *= 2.0
            rejected(model, "YOUNG_MODULUS")

            model, part, _, element = helper._build(geometry, law)
            part.CreateNewNode(5, 2.0, 2.0, 0.0)
            rejected(model, "node counts")

            model, part, _, element = helper._build(geometry, law)
            part.CreateNewElement(
                geometry["element"], 2, [1, 2, 3, 4], part.Properties[1])
            rejected(model, "element counts")

            temperature_file = tmp / "temperature_only.json"
            self._save(export_model, temperature_file, False, True, False)
            model, part, _, _ = helper._build(geometry, law)
            with self.assertRaisesRegex(RuntimeError, "does not contain initial_stresses"):
                self._load(model, temperature_file, True, False, False)

            fresh_model, fresh_part, _, _ = helper._build(geometry, law)
            with self.assertRaisesRegex(RuntimeError, "no stored InitialState"):
                self._save(fresh_model, tmp / "missing_stored.json", True, False, False)

            corrupted = json.loads(path.read_text())
            corrupted["version"] = 999
            corrupted_path = tmp / "corrupted.json"
            corrupted_path.write_text(json.dumps(corrupted))
            with self.assertRaisesRegex(RuntimeError, "Unsupported"):
                dam_initial_state_file.ReadFile(str(corrupted_path))

            corrupted = json.loads(path.read_text())
            first_key = next(iter(corrupted["initial_stresses"]))
            corrupted["initial_stresses"][first_key][0] = [1.0, 2.0]
            corrupted_path.write_text(json.dumps(corrupted))
            with self.assertRaisesRegex(RuntimeError, "stress-vector size"):
                dam_initial_state_file.ReadFile(str(corrupted_path))

            corrupted = json.loads(temperature_file.read_text())
            first_key = next(iter(corrupted["temperature"]))
            corrupted["temperature"]["not-an-id"] = corrupted["temperature"].pop(
                first_key)
            corrupted_path = tmp / "corrupted_key.json"
            corrupted_path.write_text(json.dumps(corrupted))
            with self.assertRaisesRegex(RuntimeError, "non-integer"):
                dam_initial_state_file.ReadFile(str(corrupted_path))

            model, part, _, _ = helper._build(geometry, law)
            self._load(model, path, True, False, False)
            with self.assertRaisesRegex(RuntimeError, "already has an initial state"):
                self._load(model, path, True, False, False)

    def test_complete_gravity_thermal_workflow(self):
        workflow = workflow_tests.TestDamSelectiveInitialization()
        with tempfile.TemporaryDirectory(prefix="dam_state_workflow_") as tmp:
            root = Path(tmp)
            stage_a = root / "stage_a.json"
            stage_b = root / "stage_b.json"
            stage_c_stored = root / "stage_c_stored.json"

            settings_a = workflow._parameters(thermal=True, increment="none")
            settings_a["problem_data"]["end_time"] = 1.0
            settings_a["loads_process_list"].append({
                "python_module": "save_dam_initial_state_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "save_initial_stresses": True,
                    "save_initial_temperature": False,
                    "save_reference_temperature": False,
                    "stress_source": "current_cauchy",
                    "output_file_name": str(stage_a)}})
            workflow._run(root / "A", settings_a, thermal=True)

            settings_b = workflow._parameters(thermal=True, increment="thermal")
            settings_b["loads_process_list"].append({
                "python_module": "load_dam_initial_state_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "input_file_name": str(stage_a),
                    "transfer_initial_stresses": True,
                    "transfer_initial_temperature": False,
                    "transfer_reference_temperature": False,
                    "reference_temperature_policy": "preserve"}})
            settings_b["loads_process_list"].append({
                "python_module": "save_dam_initial_state_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "save_initial_stresses": True,
                    "save_initial_temperature": True,
                    "save_reference_temperature": True,
                    "stress_source": "stored_initial_stress",
                    "output_file_name": str(stage_b)}})
            result_b = workflow._run(root / "B", settings_b, thermal=True)

            settings_c = workflow._parameters(
                thermal=True, increment="horizontal", temperature_offset=5.0)
            settings_c["constraints_process_list"] = settings_c[
                "constraints_process_list"][1:]
            settings_c["loads_process_list"].append({
                "python_module": "load_dam_initial_state_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "input_file_name": str(stage_b),
                    "transfer_initial_stresses": True,
                    "transfer_initial_temperature": True,
                    "transfer_reference_temperature": False,
                    "reference_temperature_policy": "initial_temperature"}})
            settings_c["loads_process_list"].append({
                "python_module": "save_dam_initial_state_process",
                "kratos_module": "KratosMultiphysics.DamApplication",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "save_initial_stresses": True,
                    "save_initial_temperature": False,
                    "save_reference_temperature": False,
                    "stress_source": "stored_initial_stress",
                    "output_file_name": str(stage_c_stored)}})
            result_c = workflow._run(root / "C", settings_c, thermal=True)

            data_a = json.loads(stage_a.read_text())
            data_b = json.loads(stage_b.read_text())
            data_c = json.loads(stage_c_stored.read_text())
            self.assertEqual(data_b["fields"],
                             ["initial_stresses", "temperature", "reference_temperature"])
            self.assertEqual(
                self._stress_map_error(data_b["initial_stresses"],
                                       data_a["initial_stresses"]), 0.0)
            self.assertEqual(
                self._stress_map_error(data_c["initial_stresses"],
                                       data_a["initial_stresses"]), 0.0)

            initial = result_c["initial"]
            self.assertEqual(
                max(abs(value) for displacement in initial["displacement"]
                    for value in displacement), 0.0)
            self.assertLess(max(abs(value) for value in initial["residual"][4:]), 1.e-10)
            self.assertEqual(
                self._max_error(initial["stress"], self._flatten_stress(data_a)),
                0.0)
            self.assertAlmostEqual(
                max(abs(value) for value in initial["residual"][4:]), 0.0, delta=1.e-10)
            for index, temperature in enumerate(result_b["steps"][-1]["temperature"]):
                self.assertEqual(initial["temperature"][index], temperature)
                self.assertEqual(initial["reference_temperature"][index], temperature)
            for step in result_c["steps"]:
                for actual, expected in zip(
                        step["reference_temperature"], initial["reference_temperature"]):
                    self.assertEqual(actual, expected)

            first_displacement = max(
                abs(value) for displacement in result_c["steps"][0]["displacement"]
                for value in displacement)
            final_displacement = max(
                abs(value) for displacement in result_c["steps"][-1]["displacement"]
                for value in displacement)
            self.assertGreater(first_displacement, 1.e-8)
            self.assertGreater(final_displacement, 10.0 * first_displacement)
            first_stress_change = max(
                abs(a - b) for a, b in zip(
                    initial["stress"][0], result_c["steps"][0]["stress"][0]))
            final_stress_change = max(
                abs(a - b) for a, b in zip(
                    initial["stress"][0], result_c["steps"][-1]["stress"][0]))
            self.assertGreater(first_stress_change, 1.e-4)
            self.assertGreater(final_stress_change, 10.0 * first_stress_change)
            self.assertAlmostEqual(
                sum(reaction[1] for reaction in result_c["steps"][-1]["reaction"]),
                40.0, delta=1.e-9)

    def test_complete_workflow_through_dam_analysis_run(self):
        """Three real DamAnalysis.Run() executions plus policy/compatibility cases."""
        workflow = workflow_tests.TestDamSelectiveInitialization()
        with tempfile.TemporaryDirectory(prefix="dam_state_full_run_") as tmp:
            root = Path(tmp)
            mechanical_file = root / "mechanical.json"
            stage_a_file = root / "stage_a.json"
            stage_b_file = root / "stage_b.json"

            # Purely mechanical gravity reference.
            settings = workflow._parameters(thermal=False, increment="none")
            settings["problem_data"]["end_time"] = 5.0
            self._add_save(settings, mechanical_file, "current_cauchy")
            self._run_full(
                root / "mechanical", settings, thermal=False, probe_stored=False)

            # Stage A: thermoelastic law, T == T_ref, gravity only.
            settings = workflow._parameters(thermal=True, increment="none")
            settings["problem_data"]["end_time"] = 5.0
            self._add_save(settings, stage_a_file, "current_cauchy")
            stage_a_probe = self._run_full(
                root / "stage_a", settings, thermal=True, probe_stored=False)

            data_mechanical = json.loads(mechanical_file.read_text())
            data_a = json.loads(stage_a_file.read_text())

            # Same gravity state with the mechanical and thermoelastic laws.
            self.assertLess(self._stress_map_error(
                data_mechanical["initial_stresses"],
                data_a["initial_stresses"]), 1.e-12)
            self.assertNotEqual(
                data_mechanical["materials"][0]["law"],
                data_a["materials"][0]["law"])

            # Stage A exports the final gravity state and nothing else.
            self.assertEqual(data_a["fields"], ["initial_stresses"])
            self.assertEqual(data_a["stress_source"], "current_cauchy")
            self.assertNotIn("temperature", data_a)
            self.assertNotIn("reference_temperature", data_a)
            self.assertLess(self._stress_map_error(
                data_a["initial_stresses"],
                stage_a_probe["steps"][-1]["cauchy_stress"]), 1.e-12)
            self.assertLess(
                stage_a_probe["steps"][-1]["free_residual_max"], 1.e-10)
            self.assertAlmostEqual(
                stage_a_probe["steps"][-1]["time"], 5.0, delta=1.e-12)

            # Stage B: import gravity prestress and apply a transient thermal load.
            settings = workflow._parameters(thermal=True, increment="thermal")
            self._add_load(settings, stage_a_file, True, False, False, "preserve")
            self._add_save(
                settings, stage_b_file, "stored_initial_stress", True, True, True)
            stage_b_probe = self._run_full(
                root / "stage_b", settings, thermal=True)

            data_b = json.loads(stage_b_file.read_text())
            initial_b = stage_b_probe["initial"]
            final_b = stage_b_probe["steps"][-1]
            self.assertEqual(initial_b["time"], 0.0)
            self.assertEqual(initial_b["step"], 0)
            self.assertEqual(
                max(abs(value) for displacement in initial_b["displacement"]
                    for value in displacement), 0.0)
            self.assertLess(initial_b["free_residual_max"], 1.e-10)
            self.assertLess(self._stress_map_error(
                initial_b["cauchy_stress"], data_a["initial_stresses"]), 1.e-12)
            self.assertLess(max(
                abs(value) for values in initial_b["thermal_strain"].values()
                for point in values for value in point), 1.e-12)
            self.assertGreater(max(
                abs(a - b) for a, b in zip(
                    initial_b["temperature"], final_b["temperature"])), 1.0)
            self.assertGreater(self._stress_map_error(
                final_b["cauchy_stress"], data_a["initial_stresses"]), 0.1)
            self.assertEqual(self._stress_map_error(
                final_b["stored_initial_stress"],
                data_a["initial_stresses"]), 0.0)
            self.assertEqual(self._stress_map_error(
                final_b["stored_initial_stress"],
                initial_b["stored_initial_stress"]), 0.0)

            self.assertEqual(data_b["fields"], [
                "initial_stresses", "temperature", "reference_temperature"])
            self.assertEqual(data_b["stress_source"], "stored_initial_stress")
            self.assertEqual(self._stress_map_error(
                data_b["initial_stresses"], data_a["initial_stresses"]), 0.0)
            self.assertGreater(self._stress_map_error(
                data_b["initial_stresses"], final_b["cauchy_stress"]), 0.1)
            for index, node_id in enumerate(final_b["nodes"]):
                self.assertEqual(
                    data_b["temperature"][str(node_id)],
                    final_b["temperature"][index])
                self.assertEqual(
                    data_b["reference_temperature"][str(node_id)],
                    final_b["reference_temperature"][index])

            # Stage C: independent transient analysis importing all fields.
            settings = workflow._parameters(
                thermal=True, increment="thermal", temperature_offset=5.0)
            settings["constraints_process_list"] = \
                settings["constraints_process_list"][1:]
            settings["loads_process_list"].append({
                "python_module": "assign_vector_variable_process",
                "kratos_module": "KratosMultiphysics",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "variable_name": "VOLUME_ACCELERATION",
                    "constrained": False,
                    "interval": [2.5, "End"],
                    "value": [2.0, -10.0, 0.0]}})
            self._add_load(
                settings, stage_b_file, True, True, False, "initial_temperature")
            stage_c_probe = self._run_full(
                root / "stage_c", settings, thermal=True)

            initial_c = stage_c_probe["initial"]
            first_c = stage_c_probe["steps"][0]
            final_c = stage_c_probe["steps"][-1]
            self.assertEqual(initial_c["time"], 0.0)
            self.assertEqual(initial_c["step"], 0)
            self.assertEqual(
                max(abs(value) for displacement in initial_c["displacement"]
                    for value in displacement), 0.0)
            self.assertAlmostEqual(
                sum(reaction[1] for reaction in initial_c["reaction"]), 0.0,
                delta=1.e-12)
            self.assertLess(initial_c["free_residual_max"], 1.e-10)
            self.assertEqual(self._stress_map_error(
                initial_c["stored_initial_stress"],
                data_a["initial_stresses"]), 0.0)
            self.assertEqual(self._stress_map_error(
                initial_c["cauchy_stress"], data_a["initial_stresses"]), 0.0)
            self.assertEqual(initial_c["temperature"], final_b["temperature"])
            self.assertEqual(
                initial_c["reference_temperature"], initial_c["temperature"])
            self.assertLess(max(
                abs(value) for values in initial_c["thermal_strain"].values()
                for point in values for value in point), 1.e-12)
            for index, (temperature_buffer, reference_buffer) in enumerate(zip(
                    initial_c["temperature_buffer"],
                    initial_c["reference_temperature_buffer"])):
                self.assertGreaterEqual(len(temperature_buffer), 2)
                for value in temperature_buffer:
                    self.assertEqual(value, initial_c["temperature"][index])
                for value in reference_buffer:
                    self.assertEqual(
                        value, initial_c["reference_temperature"][index])
            for step in stage_c_probe["steps"]:
                self.assertEqual(
                    step["reference_temperature"],
                    initial_c["reference_temperature"])
            for index, buffer in enumerate(first_c["temperature_buffer"]):
                self.assertEqual(buffer[1], initial_c["temperature"][index])
                self.assertTrue(all(math.isfinite(value) for value in buffer))
            self.assertGreater(max(
                abs(buffer[0] - buffer[1])
                for buffer in first_c["temperature_buffer"]), 1.e-6)

            first_displacement = max(
                abs(value) for displacement in first_c["displacement"]
                for value in displacement)
            final_displacement = max(
                abs(value) for displacement in final_c["displacement"]
                for value in displacement)
            self.assertGreater(first_displacement, 1.e-8)
            self.assertGreater(final_displacement, 10.0 * first_displacement)
            self.assertGreater(max(
                abs(a - b) for a, b in zip(
                    final_c["temperature"], initial_c["temperature"])), 1.0)
            self.assertGreater(self._stress_map_error(
                final_c["cauchy_stress"], data_a["initial_stresses"]), 1.0)
            self.assertEqual(self._stress_map_error(
                final_c["stored_initial_stress"],
                data_a["initial_stresses"]), 0.0)
            self.assertAlmostEqual(
                sum(reaction[1] for reaction in final_c["reaction"]), 40.0,
                delta=1.e-8)

            # preserve: imported temperature fields win over the existing
            # temperature/reference processes and T_ref remains constant.
            preserve_file = root / "stage_b_preserve.json"
            preserve_data = json.loads(stage_b_file.read_text())
            preserve_data["reference_temperature"] = {
                key: value + 10.0
                for key, value in preserve_data["reference_temperature"].items()}
            preserve_file.write_text(json.dumps(preserve_data))
            settings = workflow._parameters(thermal=True, increment="none")
            settings["constraints_process_list"] = \
                settings["constraints_process_list"][1:]
            self._add_load(settings, preserve_file, False, True, True, "preserve")
            preserve_probe = self._run_full(
                root / "preserve", settings, thermal=True, probe_stored=False)
            preserved_initial = preserve_probe["initial"]
            for index, node_id in enumerate(preserved_initial["nodes"]):
                self.assertEqual(
                    preserved_initial["temperature"][index],
                    data_b["temperature"][str(node_id)])
                self.assertEqual(
                    preserved_initial["reference_temperature"][index],
                    preserve_data["reference_temperature"][str(node_id)])
            for step in preserve_probe["steps"]:
                self.assertEqual(
                    step["reference_temperature"],
                    preserved_initial["reference_temperature"])

            # An active reference-temperature process conflicts with importing
            # the reference field and is rejected instead of silently winning.
            settings = workflow._parameters(thermal=True, increment="none")
            self._add_load(settings, preserve_file, False, True, True, "preserve")
            self._run_full(
                root / "preserve_conflict", settings, thermal=True,
                expected_error="Cannot import NODAL_REFERENCE_TEMPERATURE")

            # preserve with stresses now supports T != T_ref. The imported
            # stored prestress is reproduced at zero displacement while both
            # temperature fields remain unchanged.
            settings = workflow._parameters(thermal=True, increment="thermal")
            settings["constraints_process_list"] = \
                settings["constraints_process_list"][1:]
            self._add_load(
                settings, stage_b_file, True, True, True, "preserve")
            preserve_stress_probe = self._run_full(
                root / "preserve_stress", settings, thermal=True)
            preserved_stress_initial = preserve_stress_probe["initial"]
            for index, node_id in enumerate(preserved_stress_initial["nodes"]):
                self.assertEqual(
                    preserved_stress_initial["temperature"][index],
                    data_b["temperature"][str(node_id)])
                self.assertEqual(
                    preserved_stress_initial["reference_temperature"][index],
                    data_b["reference_temperature"][str(node_id)])
            self.assertEqual(
                max(abs(value) for displacement in
                    preserved_stress_initial["displacement"]
                    for value in displacement), 0.0)
            self.assertEqual(self._stress_map_error(
                preserved_stress_initial["cauchy_stress"],
                data_b["initial_stresses"]), 0.0)
            self.assertGreater(max(
                abs(value) for values in
                preserved_stress_initial["thermal_strain"].values()
                for point in values for value in point), 1.e-8)

            # Mechanical -> thermoelastic transfer is rejected, not relaxed.
            settings = workflow._parameters(thermal=True, increment="thermal")
            self._add_load(
                settings, mechanical_file, True, False, False, "preserve")
            self._run_full(
                root / "mismatch", settings, thermal=True,
                expected_error="constitutive law mismatch")

    def test_total_stress_import_with_independent_temperatures(self):
        """A->B->C total-stress transfer with T != T_ref and preserve policy."""
        workflow = workflow_tests.TestDamSelectiveInitialization()
        with tempfile.TemporaryDirectory(prefix="dam_total_stress_") as tmp:
            root = Path(tmp)
            stage_a_file = root / "stage_a.json"
            stage_b_file = root / "stage_b.json"

            # Stage A: nonuniform T and T_ref, T != T_ref, gravity plus thermal
            # stress. Export the equilibrated total Cauchy stress.
            settings_a = workflow._parameters(
                thermal=True, increment="none", temperature_offset=5.0)
            settings_a["problem_data"]["end_time"] = 1.0
            self._add_save(
                settings_a, stage_a_file, "current_cauchy", True, True, True)
            stage_a_probe = self._run_full(
                root / "A", settings_a, thermal=True, probe_stored=False,
                temperature_offset=5.0)

            data_a = json.loads(stage_a_file.read_text())
            final_a = stage_a_probe["steps"][-1]
            self.assertEqual(data_a["fields"], [
                "initial_stresses", "temperature", "reference_temperature"])
            self.assertEqual(data_a["stress_source"], "current_cauchy")
            self.assertAlmostEqual(final_a["time"], 1.0, delta=1.e-12)
            self.assertLess(final_a["free_residual_max"], 1.e-10)
            self.assertGreater(max(
                abs(temperature - reference) for temperature, reference in zip(
                    final_a["temperature"], final_a["reference_temperature"])),
                1.0)
            self.assertEqual(self._stress_map_error(
                data_a["initial_stresses"], final_a["cauchy_stress"]), 0.0)
            for index, node_id in enumerate(final_a["nodes"]):
                self.assertEqual(
                    data_a["temperature"][str(node_id)],
                    final_a["temperature"][index])
                self.assertEqual(
                    data_a["reference_temperature"][str(node_id)],
                    final_a["reference_temperature"][index])

            # Stage B: preserve both imported temperature fields, compensate the
            # initial thermal strain, then apply a +5 K boundary increment.
            settings_b = workflow._parameters(
                thermal=True, increment="none", temperature_offset=5.0)
            settings_b["constraints_process_list"] = \
                settings_b["constraints_process_list"][1:]
            base = "20+2*y+5"
            settings_b["constraints_process_list"][0]["Parameters"]["value"] = (
                base + "+5 if t>0.5 else " + base)
            settings_b["problem_data"]["end_time"] = 3.0
            self._add_load(
                settings_b, stage_a_file, True, True, True, "preserve")
            self._add_save(
                settings_b, stage_b_file, "current_cauchy", True, True, True)
            stage_b_probe = self._run_full(
                root / "B", settings_b, thermal=True,
                temperature_offset=5.0)

            data_b = json.loads(stage_b_file.read_text())
            initial_b = stage_b_probe["initial"]
            final_b = stage_b_probe["steps"][-1]
            self.assertEqual(initial_b["temperature"], final_a["temperature"])
            self.assertEqual(
                initial_b["reference_temperature"],
                final_a["reference_temperature"])
            self.assertGreater(max(
                abs(temperature - reference) for temperature, reference in zip(
                    initial_b["temperature"],
                    initial_b["reference_temperature"])),
                1.0)
            self.assertEqual(
                max(abs(value) for displacement in initial_b["displacement"]
                    for value in displacement), 0.0)
            self.assertLess(initial_b["free_residual_max"], 1.e-10)
            self.assertEqual(self._stress_map_error(
                initial_b["cauchy_stress"], final_a["cauchy_stress"]), 0.0)
            self.assertGreater(max(
                abs(value) for values in initial_b["thermal_strain"].values()
                for point in values for value in point), 1.e-8)
            self.assertEqual(
                final_b["reference_temperature"],
                initial_b["reference_temperature"])
            self.assertGreater(max(
                abs(a - b) for a, b in zip(
                    final_b["temperature"], initial_b["temperature"])), 1.0)
            self.assertGreater(self._stress_map_error(
                final_b["cauchy_stress"], initial_b["cauchy_stress"]), 0.1)
            self.assertEqual(data_b["stress_source"], "current_cauchy")
            for index, node_id in enumerate(final_b["nodes"]):
                self.assertEqual(
                    data_b["temperature"][str(node_id)],
                    final_b["temperature"][index])
                self.assertEqual(
                    data_b["reference_temperature"][str(node_id)],
                    final_b["reference_temperature"][index])

            # Continuous reference: identical source at t=1, then the same +5 K
            # increment. Compare increments after the transfer instant.
            settings_r = workflow._parameters(
                thermal=True, increment="none", temperature_offset=5.0)
            base = "20+2*y+5"
            settings_r["constraints_process_list"][1]["Parameters"]["value"] = (
                base + "+5 if t>1.5 else " + base)
            settings_r["problem_data"]["end_time"] = 4.0
            reference_probe = self._run_full(
                root / "Reference", settings_r, thermal=True,
                probe_stored=False, temperature_offset=5.0)
            reference_at_one = next(
                step for step in reference_probe["steps"]
                if abs(step["time"] - 1.0) < 1.e-12)
            self.assertEqual(
                reference_at_one["temperature"], final_a["temperature"])
            self.assertEqual(
                reference_at_one["reference_temperature"],
                final_a["reference_temperature"])
            self.assertEqual(self._stress_map_error(
                reference_at_one["cauchy_stress"], final_a["cauchy_stress"]), 0.0)
            self.assertEqual(
                self._max_error(
                    reference_at_one["displacement"], final_a["displacement"]),
                0.0)

            def stress_increment(values, base_values):
                return {
                    key: [[a - b for a, b in zip(point, base_point)]
                          for point, base_point in zip(
                              values[key], base_values[key])]
                    for key in values}

            for step in stage_b_probe["steps"]:
                reference_step = next(
                    reference for reference in reference_probe["steps"]
                    if abs(reference["time"] - (step["time"] + 1.0)) < 1.e-12)
                displacement_increment = [
                    [a - b for a, b in zip(displacement, base_displacement)]
                    for displacement, base_displacement in zip(
                        step["displacement"], initial_b["displacement"])]
                reference_displacement_increment = [
                    [a - b for a, b in zip(displacement, base_displacement)]
                    for displacement, base_displacement in zip(
                        reference_step["displacement"],
                        reference_at_one["displacement"])]
                self.assertLess(self._max_error(
                    displacement_increment,
                    reference_displacement_increment), 1.e-10)
                self.assertLess(self._stress_map_error(
                    stress_increment(
                        step["cauchy_stress"], initial_b["cauchy_stress"]),
                    stress_increment(
                        reference_step["cauchy_stress"],
                        reference_at_one["cauchy_stress"])), 1.e-10)

            # Stage C: import the Stage B total stress and both temperature
            # fields, then apply a mechanical increment.
            settings_c = workflow._parameters(
                thermal=True, increment="none", temperature_offset=5.0)
            settings_c["constraints_process_list"] = \
                settings_c["constraints_process_list"][1:]
            settings_c["constraints_process_list"][0]["Parameters"]["value"] = (
                "20+2*y+5+5")
            settings_c["problem_data"]["end_time"] = 2.0
            settings_c["loads_process_list"].append({
                "python_module": "assign_vector_variable_process",
                "kratos_module": "KratosMultiphysics",
                "Parameters": {
                    "model_part_name": "MainModelPart.Parts_Concrete",
                    "variable_name": "VOLUME_ACCELERATION",
                    "constrained": False,
                    "interval": [0.5, "End"],
                    "value": [2.0, -10.0, 0.0]}})
            self._add_load(
                settings_c, stage_b_file, True, True, True, "preserve")
            stage_c_probe = self._run_full(
                root / "C", settings_c, thermal=True,
                temperature_offset=5.0)

            initial_c = stage_c_probe["initial"]
            final_c = stage_c_probe["steps"][-1]
            self.assertEqual(initial_c["temperature"], final_b["temperature"])
            self.assertEqual(
                initial_c["reference_temperature"],
                final_b["reference_temperature"])
            self.assertGreater(max(
                abs(temperature - reference) for temperature, reference in zip(
                    initial_c["temperature"],
                    initial_c["reference_temperature"])),
                1.0)
            self.assertEqual(
                max(abs(value) for displacement in initial_c["displacement"]
                    for value in displacement), 0.0)
            self.assertLess(initial_c["free_residual_max"], 1.e-10)
            self.assertEqual(self._stress_map_error(
                initial_c["cauchy_stress"], data_b["initial_stresses"]), 0.0)
            self.assertEqual(self._stress_map_error(
                initial_c["stored_initial_stress"],
                data_b["initial_stresses"]), 0.0)
            self.assertGreater(max(
                abs(value) for values in initial_c["thermal_strain"].values()
                for point in values for value in point), 1.e-8)
            self.assertEqual(
                final_c["reference_temperature"],
                initial_c["reference_temperature"])
            self.assertGreater(
                max(abs(value) for displacement in final_c["displacement"]
                    for value in displacement), 1.e-8)
            self.assertGreater(self._stress_map_error(
                final_c["cauchy_stress"], initial_c["cauchy_stress"]), 1.e-4)
            self.assertEqual(self._stress_map_error(
                final_c["stored_initial_stress"],
                data_b["initial_stresses"]), 0.0)
            self.assertAlmostEqual(
                sum(reaction[1] for reaction in final_c["reaction"]), 40.0,
                delta=1.e-8)

    @staticmethod
    def _flatten_stress(data):
        result = []
        for element_id in sorted(data["initial_stresses"], key=int):
            result.extend(data["initial_stresses"][element_id])
        return result


if __name__ == "__main__":
    KratosUnittest.main()
