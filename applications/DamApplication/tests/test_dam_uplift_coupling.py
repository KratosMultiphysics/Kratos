import math
import os

import KratosMultiphysics as KM
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

_FIXTURE_DIR = "uplift_coupling"
_BASE_NODE_IDS = (1, 2, 3, 4, 5)
_FIXED_NODE_IDS = (101, 102, 103, 104, 105)
_JOINT_ELEMENT_IDS = (201, 202, 203, 204)

# Reference loads of the fixture.
_TOP_LOAD = 1000.0      # per unit length, downward
_WATER_LOAD = 1000.0    # per unit length, downstream
_WATER_LEVEL = 0.2      # specific weight 10000 -> p0 = 2000 Pa
_SPECIFIC_WEIGHT = 10000.0
_LENGTH = 4.0


class _RecordingDamAnalysis(DamAnalysis):
    """DamAnalysis that keeps the converged coupled state of every step."""

    def __init__(self, model, project_parameters):
        self.coupled_history = []
        super().__init__(model, project_parameters)

    def _SolveWithUpliftCoupling(self):
        super()._SolveWithUpliftCoupling()
        self.coupled_history.append(self._CoupledState())

    def _CoupledState(self):
        model_part = self.main_model_part
        base_nodes = sorted(
            [node for node in model_part.Nodes if node.Id in _BASE_NODE_IDS],
            key=lambda node: node.Id)
        fixed_nodes = sorted(
            [node for node in model_part.Nodes if node.Id in _FIXED_NODE_IDS],
            key=lambda node: node.Id)
        openings = []
        for element_id in _JOINT_ELEMENT_IDS:
            openings.extend(model_part.GetElement(element_id).CalculateOnIntegrationPoints(
                KratosPoro.JOINT_WIDTH, model_part.ProcessInfo))

        process = (self.uplift_coupling_data[0]["process"]
                   if self.uplift_coupling_data else None)
        return {
            "time": model_part.ProcessInfo[KM.TIME],
            "pressures": [node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE)
                          for node in base_nodes],
            "openings": openings,
            "x": [node.X for node in base_nodes],
            "front": process.GetJointPosition() if process is not None else None,
            "iterations": self.GetLastUpliftCouplingIterations(),
            "reaction_y": sum(
                node.GetSolutionStepValue(KM.REACTION)[1] for node in fixed_nodes),
        }


class _ThermalSolverProxy:
    """Counts thermal solves while delegating everything to the real strategy."""

    def __init__(self, strategy, owner):
        self._strategy = strategy
        self._owner = owner

    def Solve(self):
        self._owner.thermal_solve_count += 1
        return self._strategy.Solve()

    def __getattr__(self, name):
        return getattr(self._strategy, name)


class _ThermalCountingAnalysis(_RecordingDamAnalysis):
    """Thermo-mechanical analysis that counts thermal solves."""

    def Initialize(self):
        super().Initialize()
        self.thermal_solve_count = 0
        self.solver.Thermal_Solver = _ThermalSolverProxy(
            self.solver.Thermal_Solver, self)


class TestDamUpliftCoupling(KratosUnittest.TestCase):
    """Hydro-mechanical coupling between joint opening and uplift pressure.

    The fixture is a 2D dam block resting on a horizontal elastic joint. The
    top load compresses the joint, the upstream line load bends the dam and the
    uplift pressure acts on the joint. The pressure depends on the open length,
    which in turn depends on the pressure: the coupled fixed point must be
    solved iteratively. A one-way (single) update gives a different solution.
    """

    # ------------------------------------------------------------------
    # Fixture helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _LoadParameters(file_name="ProjectParameters.json"):
        fixture = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), _FIXTURE_DIR,
            file_name)
        with open(fixture) as parameter_file:
            return KM.Parameters(parameter_file.read())

    @staticmethod
    def _ConditionLoadEntry(parameters, part_name):
        for i in range(parameters["loads_process_list"].size()):
            process = parameters["loads_process_list"][i]
            if part_name in process["Parameters"]["model_part_name"].GetString():
                return process["Parameters"]
        raise RuntimeError("No load process found for '{}'".format(part_name))

    def _SetConditionLoad(self, parameters, part_name, modulus):
        self._ConditionLoadEntry(parameters, part_name)["modulus"].SetDouble(modulus)

    def _SetConditionLoadExpression(self, parameters, part_name, expression):
        values = self._ConditionLoadEntry(parameters, part_name)
        values.RemoveValue("modulus")
        values.AddEmptyValue("modulus").SetString(expression)

    def _SetWaterLevel(self, parameters, level):
        for i in range(parameters["loads_process_list"].size()):
            process = parameters["loads_process_list"][i]
            if "Uplift" in process["process_name"].GetString():
                process["Parameters"]["Water_level"].SetDouble(level)
                return
        raise RuntimeError("No uplift process found")

    def _SetCouplingSetting(self, parameters, name, value):
        for i in range(parameters["loads_process_list"].size()):
            process = parameters["loads_process_list"][i]
            if "Uplift" not in process["process_name"].GetString():
                continue
            coupling = process["Parameters"]["hydro_mechanical_coupling"]
            if not coupling.Has(name):
                coupling.AddEmptyValue(name)
            if isinstance(value, bool):
                coupling[name].SetBool(value)
            elif isinstance(value, int):
                coupling[name].SetInt(value)
            else:
                coupling[name].SetDouble(value)
            return
        raise RuntimeError("No uplift process found")

    def _Run(self, parameters, analysis_class=DamAnalysis):
        model = KM.Model()
        analysis = analysis_class(model, parameters)
        analysis.Run()
        return analysis

    def _BaselineParameters(self):
        parameters = self._LoadParameters()
        self._SetConditionLoad(parameters, "TopLoad", _TOP_LOAD)
        self._SetConditionLoad(parameters, "WaterLoad", _WATER_LOAD)
        self._SetWaterLevel(parameters, _WATER_LEVEL)
        return parameters

    # ------------------------------------------------------------------
    # Result helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _Resultant(state):
        return sum(
            0.5 * (state["x"][i + 1] - state["x"][i]) *
            (state["pressures"][i + 1] + state["pressures"][i])
            for i in range(len(state["x"]) - 1))

    def _AssertAnalysisConsistent(self, analysis, expected_front,
                                  pressure_delta=1.0e-10):
        """The stored pressure must be the one produced by the final opening."""
        process = analysis.uplift_coupling_data[0]["process"]
        model_part = analysis.main_model_part
        before = {node.Id: node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE)
                  for node in model_part.Nodes if node.Id in _BASE_NODE_IDS}

        # Re-evaluate the process from the current (converged) opening.
        process.ExecuteInitializeSolutionStep()

        after = {node.Id: node.GetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE)
                 for node in model_part.Nodes if node.Id in _BASE_NODE_IDS}
        for node_id in before:
            self.assertAlmostEqual(
                before[node_id], after[node_id], delta=pressure_delta)
        self.assertAlmostEqual(process.GetMaxJointWidthChange(), 0.0, delta=1.0e-14)
        self.assertAlmostEqual(process.GetJointPosition(), expected_front, delta=1.0e-12)

    # ------------------------------------------------------------------
    # Benchmark A: closed joint
    # ------------------------------------------------------------------

    def test_closed_joint_does_not_propagate(self):
        parameters = self._LoadParameters()
        self._SetConditionLoad(parameters, "TopLoad", 4000.0)
        self._SetConditionLoad(parameters, "WaterLoad", 0.0)
        self._SetWaterLevel(parameters, 0.05)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            analysis = self._Run(parameters, _RecordingDamAnalysis)

        state = analysis.coupled_history[-1]
        self.assertEqual(len(analysis.coupled_history), 1)
        self.assertAlmostEqual(state["front"], 0.0, delta=1.0e-12)
        self.assertEqual(state["iterations"], 1)

        # p0 = specific weight * water level = 500 Pa; uncracked linear decay.
        expected = [500.0, 375.0, 250.0, 125.0, 0.0]
        for actual, reference in zip(state["pressures"], expected):
            self.assertAlmostEqual(actual, reference, delta=1.0e-9)
        for width in state["openings"]:
            self.assertAlmostEqual(width, 0.0, delta=1.0e-15)

        uplift = self._Resultant(state)
        self.assertAlmostEqual(uplift, 0.5 * 500.0 * _LENGTH, delta=1.0e-9)
        # Vertical equilibrium: top load - uplift + reaction = 0.
        self.assertAlmostEqual(state["reaction_y"], 4000.0 * _LENGTH - uplift,
                               delta=1.0e-9)

        self._AssertAnalysisConsistent(analysis, 0.0)

    # ------------------------------------------------------------------
    # Benchmark B: coupled crack propagation (the key regression)
    # ------------------------------------------------------------------

    def test_coupled_crack_propagation_differs_from_one_way(self):
        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            coupled = self._Run(self._BaselineParameters(), _RecordingDamAnalysis)

            one_way_parameters = self._BaselineParameters()
            self._SetCouplingSetting(one_way_parameters, "enabled", False)
            one_way = self._Run(one_way_parameters, _RecordingDamAnalysis)

        coupled_state = coupled.coupled_history[-1]
        one_way_state = one_way._CoupledState()

        # The converged crack front is at x = 3.
        self.assertAlmostEqual(coupled_state["front"], 3.0, delta=1.0e-12)
        self.assertGreaterEqual(coupled_state["iterations"], 2)

        # Full hydrostatic pressure over the open length and zero at the toe.
        expected = [2000.0, 2000.0, 2000.0, 2000.0, 0.0]
        for actual, reference in zip(coupled_state["pressures"], expected):
            self.assertAlmostEqual(actual, reference, delta=1.0e-9)

        uplift = self._Resultant(coupled_state)
        self.assertAlmostEqual(uplift, 2000.0 * 3.0 + 0.5 * 2000.0 * 1.0,
                               delta=1.0e-9)
        self.assertAlmostEqual(
            coupled_state["reaction_y"], _TOP_LOAD * _LENGTH - uplift,
            delta=1.0e-9)

        # A single (one-way) update is not enough: the pressure used by the
        # mechanical solve must be the converged one.
        self.assertIsNone(one_way_state["front"])
        one_way_pressures = [2000.0, 1500.0, 1000.0, 500.0, 0.0]
        for actual, reference in zip(one_way_state["pressures"], one_way_pressures):
            self.assertAlmostEqual(actual, reference, delta=1.0e-9)
        one_way_uplift = self._Resultant(one_way_state)
        self.assertAlmostEqual(one_way_uplift, 4000.0, delta=1.0e-9)
        # Global vertical equilibrium of the one-way solution: the reaction
        # balances the top load minus the uplift actually applied.
        self.assertAlmostEqual(
            one_way_state["reaction_y"],
            _TOP_LOAD * _LENGTH - one_way_uplift, delta=1.0e-9)

        # The feedback must change the result substantially.
        pressure_differences = [
            abs(converged - one_way_value)
            for converged, one_way_value in zip(
                coupled_state["pressures"], one_way_state["pressures"])]
        # The upstream toe already carries the full pressure in both cases; the
        # feedback must change the rest of the base substantially.
        self.assertGreater(max(pressure_differences), 1000.0)
        self.assertGreater(
            max(abs(a - b) for a, b in zip(
                coupled_state["openings"], one_way_state["openings"])),
            0.4 * max(coupled_state["openings"]))
        self.assertGreater(
            max(coupled_state["openings"]),
            1.5 * max(one_way_state["openings"]))

        self._AssertAnalysisConsistent(coupled, 3.0)

    # ------------------------------------------------------------------
    # Benchmark C: convergence and tolerance stability
    # ------------------------------------------------------------------

    def test_coupling_convergence_and_tolerance_stability(self):
        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            default_parameters = self._BaselineParameters()
            default_run = self._Run(default_parameters, _RecordingDamAnalysis)

            tight_parameters = self._BaselineParameters()
            self._SetCouplingSetting(tight_parameters, "pressure_relative_tolerance", 1.0e-12)
            self._SetCouplingSetting(tight_parameters, "pressure_absolute_tolerance", 1.0e-10)
            self._SetCouplingSetting(tight_parameters, "joint_width_absolute_tolerance", 1.0e-12)
            self._SetCouplingSetting(tight_parameters, "front_relative_tolerance", 1.0e-12)
            tight_run = self._Run(tight_parameters, _RecordingDamAnalysis)

        default_state = default_run.coupled_history[-1]
        tight_state = tight_run.coupled_history[-1]

        self.assertGreaterEqual(default_state["iterations"], 2)
        self.assertGreaterEqual(tight_state["iterations"], 2)

        self.assertAlmostEqual(default_state["front"], tight_state["front"], delta=1.0e-12)
        for a, b in zip(default_state["pressures"], tight_state["pressures"]):
            self.assertAlmostEqual(a, b, delta=1.0e-6)
        for a, b in zip(default_state["openings"], tight_state["openings"]):
            self.assertAlmostEqual(a, b, delta=1.0e-9)

        self._AssertAnalysisConsistent(default_run, 3.0)
        self._AssertAnalysisConsistent(tight_run, 3.0)

    # ------------------------------------------------------------------
    # Benchmark D: opening, closure and reopening
    # ------------------------------------------------------------------

    def test_closure_and_reopening_sequence(self):
        parameters = self._BaselineParameters()
        parameters["problem_data"]["end_time"].SetDouble(3.0)
        parameters["problem_data"]["time_step"].SetDouble(1.0)
        # Low at t=1 (open), high at t=2 (closed), low at t=3 (open again).
        self._SetConditionLoadExpression(
            parameters, "TopLoad", "2500.0 + 1500.0*cos(3.141592653589793*t)")

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            analysis = self._Run(parameters, _RecordingDamAnalysis)

        history = analysis.coupled_history
        self.assertEqual(len(history), 3)

        opened, closed, reopened = history

        # 1) Opened joint.
        self.assertAlmostEqual(opened["front"], 3.0, delta=1.0e-12)
        for width in opened["openings"]:
            self.assertGreater(width, -1.0e-15)
        self.assertGreater(max(opened["openings"]), 0.0)

        # 2) Mechanically closed joint: the pressure returns to the uncracked
        #    linear distribution and the aperture is zero everywhere.
        self.assertAlmostEqual(closed["front"], 0.0, delta=1.0e-12)
        expected_closed = [2000.0, 1500.0, 1000.0, 500.0, 0.0]
        for actual, reference in zip(closed["pressures"], expected_closed):
            self.assertAlmostEqual(actual, reference, delta=1.0e-9)
        for width in closed["openings"]:
            self.assertAlmostEqual(width, 0.0, delta=1.0e-15)

        # 3) Reopening is reversible for the elastic cohesive law: the same
        #    aperture and pressure as in the first step are recovered.
        self.assertAlmostEqual(reopened["front"], 3.0, delta=1.0e-12)
        for a, b in zip(reopened["pressures"], opened["pressures"]):
            self.assertAlmostEqual(a, b, delta=1.0e-9)
        for a, b in zip(reopened["openings"], opened["openings"]):
            self.assertAlmostEqual(a, b, delta=1.0e-9)

    # ------------------------------------------------------------------
    # Benchmark E: non-convergence must be reported
    # ------------------------------------------------------------------

    def test_non_convergence_is_reported(self):
        parameters = self._BaselineParameters()
        self._SetCouplingSetting(parameters, "max_iterations", 1)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            with self.assertRaises(RuntimeError) as context:
                self._Run(parameters)

        message = str(context.exception)
        self.assertIn("did not converge", message)
        self.assertIn("NOT advanced", message)

    # ------------------------------------------------------------------
    # Benchmark F: pressure relaxation
    # ------------------------------------------------------------------

    def test_pressure_relaxation_factor(self):
        parameters = self._BaselineParameters()
        self._SetCouplingSetting(parameters, "relaxation_factor", 0.5)
        self._SetCouplingSetting(parameters, "echo_level", 2)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            analysis = self._Run(parameters, _RecordingDamAnalysis)

        history = analysis.GetUpliftCouplingPressureHistory()
        self.assertGreater(len(history), 1)

        # The first mechanical solve uses the uncracked (one-way) pressure.
        # The pressure recomputed from the resulting opening is the full
        # hydrostatic distribution, so the applied pressure is the relaxed
        # value between both.
        used, raw, applied = history[0][0]
        expected_used = {1: 2000.0, 2: 1500.0, 3: 1000.0, 4: 500.0, 5: 0.0}
        expected_raw = {1: 2000.0, 2: 2000.0, 3: 2000.0, 4: 2000.0, 5: 0.0}
        for node_id in expected_used:
            self.assertAlmostEqual(used[node_id], expected_used[node_id],
                                   delta=1.0e-9)
            self.assertAlmostEqual(raw[node_id], expected_raw[node_id],
                                   delta=1.0e-9)
            self.assertAlmostEqual(
                applied[node_id],
                0.5 * used[node_id] + 0.5 * raw[node_id], delta=1.0e-12)
        self.assertAlmostEqual(applied[2], 1750.0, delta=1.0e-12)

        # Every iteration relaxes between the pressure used by the mechanical
        # solve and the one recomputed from the resulting joint state. The
        # pressure used must be the one applied in the previous iteration.
        previous_applied = applied
        for record in history[1:]:
            used, raw, applied = record[0]
            for node_id in used:
                self.assertAlmostEqual(
                    used[node_id], previous_applied[node_id], delta=1.0e-12)
                self.assertAlmostEqual(
                    applied[node_id],
                    0.5 * used[node_id] + 0.5 * raw[node_id], delta=1.0e-12)
            previous_applied = applied

        # The last raw pressure (recomputed from the final opening) and the
        # pressure used by the last mechanical solve agree within the
        # convergence criterion.
        used, raw, applied = history[-1][0]
        settings = analysis.uplift_coupling_data[0]["settings"]
        pressure_tolerance = (
            settings["pressure_absolute_tolerance"].GetDouble() +
            settings["pressure_relative_tolerance"].GetDouble() *
            max(abs(value) for value in raw.values()))
        for node_id in used:
            self.assertLessEqual(
                abs(raw[node_id] - used[node_id]), pressure_tolerance)

        # The stored (relaxed) pressure is consistent with the final opening
        # within the convergence tolerance.
        self._AssertAnalysisConsistent(
            analysis, 3.0, pressure_delta=pressure_tolerance)

        # Relaxation changes the path, not the fixed point: the converged
        # solution is the same full hydrostatic distribution.
        state = analysis.coupled_history[-1]
        self.assertAlmostEqual(state["front"], 3.0, delta=1.0e-9)
        for actual, reference in zip(
                state["pressures"], [2000.0, 2000.0, 2000.0, 2000.0, 0.0]):
            self.assertAlmostEqual(actual, reference, delta=pressure_tolerance)

        # The final mechanical equilibrium corresponds to the applied
        # pressure. The relaxation residue is bounded by the tolerance.
        uplift = self._Resultant(state)
        self.assertAlmostEqual(
            state["reaction_y"], _TOP_LOAD * _LENGTH - uplift,
            delta=_LENGTH * pressure_tolerance)

    def test_relaxed_coupling_reports_non_convergence(self):
        parameters = self._BaselineParameters()
        self._SetCouplingSetting(parameters, "relaxation_factor", 0.5)
        self._SetCouplingSetting(parameters, "max_iterations", 2)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            with self.assertRaises(RuntimeError) as context:
                self._Run(parameters)

        message = str(context.exception)
        self.assertIn("did not converge", message)
        self.assertIn("NOT advanced", message)

    # ------------------------------------------------------------------
    # Benchmark G: thermo-mechanical path
    # ------------------------------------------------------------------

    def test_thermo_mechanical_coupling_solves_thermal_once_per_step(self):
        parameters = self._LoadParameters("ProjectParameters_thermo.json")
        parameters["problem_data"]["end_time"].SetDouble(3.0)
        parameters["problem_data"]["time_step"].SetDouble(1.0)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            analysis = self._Run(parameters, _ThermalCountingAnalysis)

        # The thermal problem does not depend on the uplift pressure: exactly
        # one thermal solve per physical time step, never one per coupling
        # iteration.
        self.assertEqual(analysis.thermal_solve_count, 3)
        self.assertEqual(len(analysis.coupled_history), 3)

        # The mechanical problem is iterated with the uplift pressure.
        self.assertGreaterEqual(analysis.coupled_history[0]["iterations"], 2)
        self.assertGreaterEqual(analysis.GetLastUpliftCouplingIterations(), 1)

        # The final coupled state matches the mechanical benchmark.
        state = analysis.coupled_history[-1]
        self.assertAlmostEqual(state["front"], 3.0, delta=1.0e-12)
        expected = [2000.0, 2000.0, 2000.0, 2000.0, 0.0]
        for actual, reference in zip(state["pressures"], expected):
            self.assertAlmostEqual(actual, reference, delta=1.0e-9)
        uplift = self._Resultant(state)
        self.assertAlmostEqual(
            state["reaction_y"], _TOP_LOAD * _LENGTH - uplift, delta=1.0e-9)
        self._AssertAnalysisConsistent(analysis, 3.0)

        # The fixture is isothermal at T = 0.
        temperatures = sorted(set(
            round(node.GetSolutionStepValue(KM.TEMPERATURE), 12)
            for node in analysis.main_model_part.Nodes))
        self.assertEqual(temperatures, [0.0])

    def test_thermo_mechanical_without_coupling_keeps_original_path(self):
        parameters = self._LoadParameters("ProjectParameters_thermo.json")
        self._SetCouplingSetting(parameters, "enabled", False)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            analysis = self._Run(parameters, _ThermalCountingAnalysis)

        self.assertFalse(analysis._use_uplift_coupling)
        self.assertEqual(analysis.thermal_solve_count, 1)
        self.assertEqual(len(analysis.coupled_history), 0)
        self.assertEqual(analysis.GetLastUpliftCouplingIterations(), 0)

        # Single update: the uncracked pressure distribution and no coupling.
        state = analysis._CoupledState()
        self.assertIsNone(state["front"])
        expected = [2000.0, 1500.0, 1000.0, 500.0, 0.0]
        for actual, reference in zip(state["pressures"], expected):
            self.assertAlmostEqual(actual, reference, delta=1.0e-9)

    # ------------------------------------------------------------------
    # Benchmark H: unsupported configurations and compatibility
    # ------------------------------------------------------------------

    def test_dynamic_analysis_with_coupling_is_rejected(self):
        parameters = self._BaselineParameters()
        parameters["solver_settings"]["mechanical_solver_settings"][
            "solution_type"].SetString("Dynamic")

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            with self.assertRaises(RuntimeError) as context:
                self._Run(parameters)

        message = str(context.exception)
        self.assertIn("Quasi-Static", message)
        self.assertIn("Dynamic", message)

    def test_arc_length_analysis_with_coupling_is_rejected(self):
        parameters = self._BaselineParameters()
        mechanical = parameters["solver_settings"]["mechanical_solver_settings"]
        mechanical["strategy_type"].SetString("Arc-Length")
        # Provide a complete arc-length configuration. The coupling guard is
        # expected to reject the analysis before the strategy is constructed.
        mechanical.AddEmptyValue("desired_iterations").SetInt(4)
        mechanical.AddEmptyValue("max_radius_factor").SetDouble(20.0)
        mechanical.AddEmptyValue("min_radius_factor").SetDouble(0.5)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            with self.assertRaises(RuntimeError) as context:
                self._Run(parameters)

        message = str(context.exception)
        self.assertIn("Newton-Raphson", message)
        self.assertIn("Arc-Length", message)

    def test_conflicting_uplift_coupling_settings_are_rejected(self):
        parameters = self._BaselineParameters()
        duplicate = None
        for i in range(parameters["loads_process_list"].size()):
            process = parameters["loads_process_list"][i]
            if "Uplift" in process["process_name"].GetString():
                duplicate = KM.Parameters(process.WriteJsonString())
                break
        self.assertIsNotNone(duplicate)
        duplicate["Parameters"]["hydro_mechanical_coupling"][
            "max_iterations"].SetInt(7)
        parameters["loads_process_list"].Append(duplicate)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            with self.assertRaises(RuntimeError) as context:
                self._Run(parameters)

        message = str(context.exception)
        self.assertIn("different", message)
        self.assertIn("hydro_mechanical_coupling", message)

    def test_analysis_without_uplift_process_uses_original_solver(self):
        parameters = self._BaselineParameters()
        kept = KM.Parameters("[]")
        for i in range(parameters["loads_process_list"].size()):
            process = parameters["loads_process_list"][i]
            if "Uplift" not in process["process_name"].GetString():
                kept.Append(process)
        parameters.RemoveValue("loads_process_list")
        parameters.AddValue("loads_process_list", kept)

        coupling_calls = []
        original = _RecordingDamAnalysis._SolveWithUpliftCoupling

        class _Probe(_RecordingDamAnalysis):
            def _SolveWithUpliftCoupling(self):
                coupling_calls.append(True)
                return original(self)

        with KratosUnittest.WorkFolderScope(_FIXTURE_DIR, __file__):
            analysis = self._Run(parameters, _Probe)

        self.assertFalse(analysis._use_uplift_coupling)
        self.assertEqual(coupling_calls, [])
        state = analysis._CoupledState()
        self.assertIsNone(state["front"])
        for pressure in state["pressures"]:
            self.assertAlmostEqual(pressure, 0.0, delta=1.0e-15)


if __name__ == "__main__":
    KratosUnittest.main()
