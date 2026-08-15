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

    def _double_call_compare(self, process, model_part, variables):
        """ExecuteBeforeSolutionLoop() twice at the same TIME and compare the
        observable state (assigned values, fixity, DOF set) after each call."""
        def snapshot():
            snap = {}
            for node in model_part.Nodes:
                nid = node.Id
                snap[nid] = {
                    "values": {v.Name(): node.GetSolutionStepValue(v) for v in variables},
                    "fixed": {v.Name(): node.IsFixed(v) for v in variables},
                    "has_dof": {v.Name(): node.HasDofFor(v) for v in variables},
                }
            return snap

        process.ExecuteBeforeSolutionLoop()
        after_first = snapshot()
        process.ExecuteBeforeSolutionLoop()
        after_second = snapshot()
        self.assertEqual(after_first, after_second,
                         "state changed between the two pre-loop calls")
        return after_first

    # ------------------------------------------------------------------
    # Idempotence of the double pre-loop callback (ExecuteBeforeSolutionLoop
    # is invoked twice around solver.Initialize in dam_analysis.py; see the
    # audit's architectural note). Each test calls it twice at the same TIME
    # and asserts the observable state is unchanged.
    # ------------------------------------------------------------------
    def test_idempotent_thermal_process(self):
        # deterministic thermal/reference process (table-driven overwrite)
        model, mp = self._make_model([KratosDam.NODAL_REFERENCE_TEMPERATURE])
        settings = self._make_wrapper_settings(
            "main.target", r"""
            {
                "model_part_name": "",
                "variable_name": "NODAL_REFERENCE_TEMPERATURE",
                "initial_value": 8.3
            }""")
        process = self._factory("impose_nodal_reference_temperature_process", settings, model)
        self._double_call_compare(process, mp, [KratosDam.NODAL_REFERENCE_TEMPERATURE])

    def test_idempotent_material_process(self):
        # material-evolution process (Young modulus from parameters, fixed)
        model, mp = self._make_model([KratosDam.NODAL_YOUNG_MODULUS])
        settings = self._make_wrapper_settings(
            "main.target", r"""
            {
                "model_part_name": "",
                "variable_name": "NODAL_YOUNG_MODULUS",
                "Young_Modulus_1": 20.0,
                "Young_Modulus_2": 30.0,
                "Young_Modulus_3": 40.0,
                "Young_Modulus_4": 50.0
            }""")
        process = self._factory("impose_nodal_young_modulus_process", settings, model)
        snap = self._double_call_compare(process, mp, [KratosDam.NODAL_YOUNG_MODULUS])
        # the wrapper documents that the scalar is automatically fixed; the DOF
        # set must be stable across both calls
        self.assertTrue(snap[1]["has_dof"]["NODAL_YOUNG_MODULUS"])

    def test_idempotent_load_process(self):
        # hydrostatic load (depth-derived pressure overwrite)
        model, mp = self._make_model(
            [KratosMultiphysics.POSITIVE_FACE_PRESSURE],
            target_submodel_name="HydroLinePressure2D_target")
        settings = self._make_wrapper_settings(
            "main.HydroLinePressure2D_target", r"""
            {
                "model_part_name": "",
                "variable_name": "POSITIVE_FACE_PRESSURE",
                "Gravity_Direction": "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction": 0.0,
                "Spe_weight": 10000.0,
                "Water_level": 10.0,
                "Water_Table": 0
            }""")
        process = self._factory("impose_water_loads_condition_process", settings, model)
        snap = self._double_call_compare(process, mp, [KratosMultiphysics.POSITIVE_FACE_PRESSURE])
        # node 1 at y=0 below water level 10 -> nonzero pressure
        self.assertGreater(snap[1]["values"]["POSITIVE_FACE_PRESSURE"], 0.0)

    def test_idempotent_random_fields(self):
        # random field: the table is generated ONCE by the wrapper (gstools);
        # the C++ process only reads the table by node id. Here the process is
        # exercised directly with a pre-built table (gstools is not required).
        model, mp = self._make_model([KratosDam.NODAL_YOUNG_MODULUS])
        table = KratosMultiphysics.PiecewiseLinearTable()
        table.AddRow(1, 3.5)
        table.AddRow(2, 7.2)
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name": "main.target",
            "variable_name": "NODAL_YOUNG_MODULUS",
            "mean_value": 0.0, "min_value": 0.0, "max_value": 10.0,
            "variance": 1.0, "corr_length": 1
        }""")
        process = KratosDam.DamRandomFieldsVariableProcess(
            mp.GetSubModelPart("target"), table, settings)
        snap = self._double_call_compare(process, mp, [KratosDam.NODAL_YOUNG_MODULUS])
        self.assertEqual(snap[1]["values"]["NODAL_YOUNG_MODULUS"], 3.5)
        self.assertEqual(snap[2]["values"]["NODAL_YOUNG_MODULUS"], 0.0)  # outside model part

    def test_idempotent_added_mass(self):
        # added-mass contribution (deterministic function of node coordinates)
        model, mp = self._make_model([KratosDam.ADDED_MASS])
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name": "main.target",
            "variable_name": "ADDED_MASS",
            "Gravity_Direction": "Y",
            "Reservoir_Bottom_Coordinate_in_Gravity_Direction": 0.0,
            "Spe_weight": 10000.0,
            "Water_level": 10.0
        }""")
        process = KratosDam.DamAddedMassConditionProcess(
            mp.GetSubModelPart("target"), settings)
        snap = self._double_call_compare(process, mp, [KratosDam.ADDED_MASS])
        self.assertGreater(snap[1]["values"]["ADDED_MASS"], 0.0)

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
    def test_nodal_young_modulus(self):
        """Positive regression for the PR #13472 nodal-Young-modulus wrapper bug.

        PR #13472 removed ``from KratosMultiphysics import *`` from
        ``impose_nodal_young_modulus_process.py`` but left the class base as the
        bare ``Process`` (now undefined), which made the module unimportable
        (``NameError``). This test verifies the one-line fix
        (``class ...(KratosMultiphysics.Process)``) end to end:
        1. module imports,
        2. factory instantiates the production wrapper,
        3. the wrapper is a valid ``KratosMultiphysics.Process``,
        4. the underlying C++ process executes on a minimal ModelPart,
        5. the expected ``NODAL_YOUNG_MODULUS`` value is assigned,
        6. fixity is exactly the documented production behaviour (the scalar is
           "automatically fixed", matching the historical ``Node::Fix``
           semantics) and no unrelated DOF/fixity changes occur.
        """
        from KratosMultiphysics.DamApplication.impose_nodal_young_modulus_process import (
            Factory as NodalYoungFactory,
        )

        model, mp = self._make_model([KratosDam.NODAL_YOUNG_MODULUS])
        settings = self._make_wrapper_settings(
            "main.target", r"""
            {
                "model_part_name": "",
                "variable_name": "NODAL_YOUNG_MODULUS",
                "Young_Modulus_1": 20.0,
                "Young_Modulus_2": 30.0,
                "Young_Modulus_3": 40.0,
                "Young_Modulus_4": 50.0
            }""")
        # the module import above already proves (1); (2) factory instantiation
        process = NodalYoungFactory(settings, model)
        self.assertIsInstance(process, KratosMultiphysics.Process)  # (3)

        # (4) the underlying process executes at the init and per-step callbacks
        process.ExecuteBeforeSolutionLoop()
        self.assertNotEqual(
            mp.GetNode(1).GetSolutionStepValue(KratosDam.NODAL_YOUNG_MODULUS), 0.0)
        mp.CloneTimeStep(1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertNotEqual(
            mp.GetNode(1).GetSolutionStepValue(KratosDam.NODAL_YOUNG_MODULUS), 0.0)

        # (5) expected value on the target node (Young_Modulus_1)
        self._assert_node_value(mp, 1, KratosDam.NODAL_YOUNG_MODULUS, 20.0)
        # node outside the process model part untouched
        self._assert_node_value(mp, 2, KratosDam.NODAL_YOUNG_MODULUS, 0.0)

        # (6) documented fixity: NODAL_YOUNG_MODULUS is fixed (wrapper comment:
        # "the scalar value is automatically fixed", Node::Fix semantics, same
        # pre/post #13472); unrelated variables must not be fixed / get DOFs.
        node = mp.GetNode(1)
        self.assertTrue(node.IsFixed(KratosDam.NODAL_YOUNG_MODULUS))
        for unrelated in (KratosMultiphysics.DISPLACEMENT_X,
                          KratosMultiphysics.DISPLACEMENT_Y,
                          KratosMultiphysics.TEMPERATURE):
            self.assertFalse(node.HasDofFor(unrelated),
                             msg="unintended DOF for %s" % unrelated.Name())
            self.assertFalse(node.IsFixed(unrelated),
                             msg="unintended fixity for %s" % unrelated.Name())

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
