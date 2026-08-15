import os
import sys
import tempfile

import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestDamProcessLifecycle(KratosUnittest.TestCase):
    """Generic process-level lifecycle harness for PR #13472 R1 processes.

    PR #13472 renamed the initialization callback of every Dam process from
    ``ExecuteInitialize`` to ``ExecuteBeforeSolutionLoop`` (coordinated with the
    analysis and the Python wrappers). The source audit classifies all of them
    L0 (lifecycle equivalent): the value is assigned before solver.Initialize
    in both versions and the consumers read it only during the solve.

    These tests lock the production-chain invariant for representative R1
    processes: through the real Python wrapper, the value is actually assigned
    at the intended initialization callback (guarding against the dead-callback
    failure mode demonstrated for Bofang) and re-applied at the per-step
    callback, touching only the target model part.
    """

    def _make_model(self, variables, target_submodel_name="target"):
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("main")
        mp.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 1.0
        for var in variables:
            mp.AddNodalSolutionStepVariable(var)
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)   # in target model part
        mp.CreateNewNode(2, 5.0, 5.0, 0.0)   # outside target model part
        sub = mp.CreateSubModelPart(target_submodel_name)
        sub.AddNodes([1])
        return model, mp
    @staticmethod
    def _factory(wrapper_module, settings, model):
        # settings = {"python_module":..., "Parameters": {...}} exactly as the
        # process_factory would hand it to the wrapper Factory
        module = __import__("KratosMultiphysics.DamApplication.%s" % wrapper_module,
                            fromlist=["Factory"])
        return module.Factory(settings, model)

    def _make_wrapper_settings(self, model_part_name, extra):
        # build {"python_module","Parameters": {...}} as the process_factory would
        params = KratosMultiphysics.Parameters(extra)
        params["model_part_name"].SetString(model_part_name)
        full = KratosMultiphysics.Parameters("""{
            "python_module": "placeholder",
            "kratos_module": "KratosMultiphysics.DamApplication",
            "process_name": "placeholder",
            "Parameters": {}
        }""")
        full["Parameters"] = params
        return full

    def _assert_node_value(self, model_part, node_id, variable, expected, tol=1.0e-9):
        self.assertAlmostEqual(model_part.GetNode(node_id).GetSolutionStepValue(variable),
                               expected, delta=tol)

    # ------------------------------------------------------------------ T
    def test_nodal_reference_temperature(self):
        model, mp = self._make_model([KratosDam.NODAL_REFERENCE_TEMPERATURE])
        settings = self._make_wrapper_settings(
            "main.target", r"""
            {
                "model_part_name": "",
                "variable_name": "NODAL_REFERENCE_TEMPERATURE",
                "initial_value": 8.3
            }""")
        process = self._factory("impose_nodal_reference_temperature_process", settings, model)

        process.ExecuteBeforeSolutionLoop()
        self._assert_node_value(mp, 1, KratosDam.NODAL_REFERENCE_TEMPERATURE, 8.3)
        self._assert_node_value(mp, 2, KratosDam.NODAL_REFERENCE_TEMPERATURE, 0.0)

        mp.CloneTimeStep(1.0)
        process.ExecuteInitializeSolutionStep()
        self._assert_node_value(mp, 1, KratosDam.NODAL_REFERENCE_TEMPERATURE, 8.3)

    def test_grouting_reference_temperature(self):
        model, mp = self._make_model([KratosDam.NODAL_REFERENCE_TEMPERATURE, KratosMultiphysics.TEMPERATURE])
        settings = self._make_wrapper_settings(
            "main.target", r"""
            {
                "model_part_name": "",
                "variable_name": "NODAL_REFERENCE_TEMPERATURE",
                "initial_value": 5.0,
                "time_grouting": 100.0
            }""")
        process = self._factory("impose_grouting_reference_temperature_process", settings, model)

        process.ExecuteBeforeSolutionLoop()
        self._assert_node_value(mp, 1, KratosDam.NODAL_REFERENCE_TEMPERATURE, 5.0)
        # node 2 is outside the process model part -> untouched
        self._assert_node_value(mp, 2, KratosDam.NODAL_REFERENCE_TEMPERATURE, 0.0)

        # at the grouting time the reference temperature is updated from TEMPERATURE
        mp.GetNode(1).SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 12.0)
        mp.CloneTimeStep(100.0)
        process.ExecuteFinalizeSolutionStep()
        self._assert_node_value(mp, 1, KratosDam.NODAL_REFERENCE_TEMPERATURE, 12.0)

    def test_reservoir_monitoring_temperature(self):
        model, mp = self._make_model(
            [KratosMultiphysics.TEMPERATURE],
            target_submodel_name="MONITORINGRESERVOIRTEMPERATURE_target")
        settings = self._make_wrapper_settings(
            "main.MONITORINGRESERVOIRTEMPERATURE_target", r"""
            {
                "model_part_name": "",
                "variable_name": "TEMPERATURE",
                "Gravity_Direction": "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction": 0.0,
                "Height_Dam": 20.0,
                "Ambient_temp": 10.0,
                "Water_level": 10.0,
                "Z_Coord_1": 0.0, "Water_temp_1": 4.0,
                "Z_Coord_2": 10.0, "Water_temp_2": 8.0,
                "Z_Coord_3": 20.0, "Water_temp_3": 12.0
            }""")
        process = self._factory("impose_reservoir_temperature_condition_process", settings, model)

        process.ExecuteBeforeSolutionLoop()
        # target node y=0, water level 10 -> assigned, below-water value
        self.assertNotEqual(mp.GetNode(1).GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0)
        # node outside the reservoir model part untouched
        self._assert_node_value(mp, 2, KratosMultiphysics.TEMPERATURE, 0.0)

        mp.CloneTimeStep(1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertNotEqual(mp.GetNode(1).GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), 0.0)

    # ------------------------------------------------------------------ M
    def test_nodal_young_modulus_wrapper_is_importable(self):
        """Canary for the OPEN PR #13472 wrapper bug.

        PR #13472 removed ``from KratosMultiphysics import *`` from
        ``impose_nodal_young_modulus_process.py`` but left the class base as the
        bare ``Process``, which is now undefined. The module therefore cannot be
        imported (``NameError``) and the process is unusable. This test pins the
        current broken behaviour (reproduction); it must be replaced by the
        functional value-assignment test once the wrapper is fixed.

        Root cause: one-line omission in #13472 — the class base was not updated
        to ``KratosMultiphysics.Process`` when the star import was removed.
        Proposed fix (reported separately, not applied here):
        ``class ImposeNodalYoungModulusProcess(KratosMultiphysics.Process):``
        """
        with self.assertRaises(NameError):
            __import__(
                "KratosMultiphysics.DamApplication.impose_nodal_young_modulus_process",
                fromlist=["Factory"])

    def test_input_table_nodal_young_modulus(self):
        model, mp = self._make_model([KratosDam.NODAL_YOUNG_MODULUS])
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("1 25.0\n")
            f.write("2 35.0\n")
            table_file = f.name
        try:
            settings = self._make_wrapper_settings(
                "main.target", r"""
                {
                    "model_part_name": "",
                    "variable_name": "NODAL_YOUNG_MODULUS",
                    "input_file_name": "",
                    "min_value": 0.0,
                    "max_value": 100.0
                }""")
            settings["Parameters"]["input_file_name"].SetString(table_file)
            process = self._factory(
                "impose_input_table_nodal_young_modulus_process", settings, model)

            process.ExecuteBeforeSolutionLoop()
            # input file present -> table values applied already at the init callback
            self._assert_node_value(mp, 1, KratosDam.NODAL_YOUNG_MODULUS, 25.0)
            # node 2 outside the process model part -> untouched
            self._assert_node_value(mp, 2, KratosDam.NODAL_YOUNG_MODULUS, 0.0)

            mp.CloneTimeStep(1.0)
            process.ExecuteInitializeSolutionStep()
            self._assert_node_value(mp, 1, KratosDam.NODAL_YOUNG_MODULUS, 25.0)
            self._assert_node_value(mp, 2, KratosDam.NODAL_YOUNG_MODULUS, 0.0)
        finally:
            os.unlink(table_file)


if __name__ == "__main__":
    KratosUnittest.main()
