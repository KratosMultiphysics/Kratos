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
    def _LoadParameters():
        fixture = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), _FIXTURE_DIR,
            "ProjectParameters.json")
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

    def _AssertAnalysisConsistent(self, analysis, expected_front):
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
            self.assertAlmostEqual(before[node_id], after[node_id], delta=1.0e-10)
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
        self.assertAlmostEqual(self._Resultant(one_way_state), 4000.0, delta=1.0e-9)

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


if __name__ == "__main__":
    KratosUnittest.main()
