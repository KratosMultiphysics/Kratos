"""Lightweight selective initial-state import/export tests."""
import json
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

            import_model, import_part, _, _ = helper._build(geometry, law)
            with self.assertRaisesRegex(RuntimeError, "T/T_ref compatibility"):
                self._load(import_model, path, True, True, True, "preserve")

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

    @staticmethod
    def _flatten_stress(data):
        result = []
        for element_id in sorted(data["initial_stresses"], key=int):
            result.extend(data["initial_stresses"][element_id])
        return result


if __name__ == "__main__":
    KratosUnittest.main()
