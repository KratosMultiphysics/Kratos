import os
import sys
import tempfile

import KratosMultiphysics
import KratosMultiphysics.DamApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

# Allow importing the experiment tooling
_EXPERIMENT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bofang_initialization_regression")
sys.path.insert(0, _EXPERIMENT_DIR)

from instrumentation.instrumented_dam_analysis import RunInstrumentedAnalysis  # noqa: E402
from comparison.analytical_bofang import BofangAnalytical  # noqa: E402


class TestBofangInitializationRegression(KratosUnittest.TestCase):
    """Regression test for the Bofang reservoir-temperature initialization.

    Invariant under test: after the full analysis initialization (and therefore
    at the beginning of the first solution step), the reservoir nodes below the
    water level must carry the analytical Bofang temperature.

    This is the physically meaningful outcome of the Bofang lifecycle. It is
    robust to the internal mechanism (ExecuteInitialize vs
    ExecuteBeforeSolutionLoop) but fails if the Bofang assignment becomes dead
    code during initialization (e.g. if the process lifecycle is reverted
    without reverting the analysis call site, PR #13472's coupling).
    """

    MONITORED_NODES = [1, 4, 7, 10, 13, 16, 17, 18]
    RESERVOIR_FACE_NODES = {1, 4, 7, 10, 13, 16}
    NODE_Y = {1: 20.0, 4: 16.0, 7: 12.0, 10: 8.0, 13: 4.0, 16: 0.0, 17: 0.0, 18: 0.0}

    def _run_and_read(self, project_file="case/ProjectParameters.json",
                      record_fixity=False):
        with KratosUnittest.WorkFolderScope("bofang_initialization_regression", __file__):
            out_dir = tempfile.mkdtemp(prefix="bofang_regression_")
            RunInstrumentedAnalysis(
                project_file,
                out_dir,
                "test",
                self.MONITORED_NODES,
                legacy_process_initialize=False,
                record_fixity=record_fixity,
            )
            rows = self._read_csv(out_dir, "lifecycle_temperatures.csv", "temperature")
            fixity = None
            if record_fixity:
                fixity = self._read_csv(out_dir, "lifecycle_fixity.csv", "is_fixed_temperature")
            return rows, fixity

    @staticmethod
    def _read_csv(out_dir, fname, value_col="temperature"):
        path = os.path.join(out_dir, fname)
        if not os.path.exists(path):
            raise FileNotFoundError("missing expected output: %s" % path)
        rows = {}
        with open(path) as f:
            header = f.readline().strip().split(",")
            for line in f:
                parts = line.strip().split(",")
                record = dict(zip(header, parts))
                stage = record["lifecycle_stage"]
                node_id = int(record["node_id"])
                rows.setdefault(stage, {})[node_id] = float(record[value_col])
        return rows

    @staticmethod
    def _check_analytical_temperature(stages, stage, node_y, reservoir_face_nodes, tol=1.0e-9):
        """Assert the analytical Bofang profile at a given lifecycle stage."""
        analytical = BofangAnalytical()
        checked = 0
        for node_id, temp in stages[stage].items():
            if node_id not in reservoir_face_nodes:
                raise AssertionError("unexpected node %d in reservoir check" % node_id)
            expected = analytical.temperature({"id": node_id, "x": 0.0, "y": node_y[node_id]})
            if expected is None:
                if abs(temp) > tol:
                    raise AssertionError("above-water node %d got a temperature %s" % (node_id, temp))
            else:
                if abs(temp - expected) > tol:
                    raise AssertionError(
                        "node %d does not match analytical Bofang: got %s expected %s"
                        % (node_id, temp, expected))
                checked += 1
        return checked

    def test_bofang_initial_temperature_matches_analytical(self):
        stages, _ = self._run_and_read()

        # stage reached right after the whole Initialize() (post-solver-init
        # ExecuteBeforeSolutionLoop), i.e. the state at the first solution step
        stage = "after_execute_before_solution_loop"
        self.assertIn(stage, stages, "lifecycle stage not recorded")

        # non-reservoir nodes must stay untouched (0.0)
        tol = 1.0e-9
        for node_id in [17, 18]:
            self.assertAlmostEqual(stages[stage][node_id], 0.0, delta=tol,
                                   msg="non-reservoir node %d got a temperature" % node_id)

        # reservoir face nodes must carry the analytical depth profile
        face = {n: stages[stage][n] for n in self.RESERVOIR_FACE_NODES}
        checked = self._check_analytical_temperature(
            {stage: face}, stage, self.NODE_Y, self.RESERVOIR_FACE_NODES)
        self.assertGreaterEqual(checked, 4)

    def test_bofang_temperature_available_before_solver_initialize(self):
        stages, _ = self._run_and_read()
        pre = stages.get("before_solver_initialize", {})
        post = stages.get("after_execute_before_solution_loop", {})
        # the water-covered face node (7: y=12 below water level 14) must carry
        # the Bofang value already before the solver is initialized
        self.assertIn(7, pre)
        self.assertAlmostEqual(pre[7], post[7], delta=1.0e-12,
                               msg="Bofang temperature changed during solver initialize")

    def test_production_wrapper_integration(self):
        """Production chain: DamAnalysis -> impose_reservoir_temperature_condition_process
        -> DamBofangConditionTemperatureProcess (C++).

        The production Python wrapper only forwards ExecuteBeforeSolutionLoop /
        ExecuteInitializeSolutionStep / ExecuteFinalizeSolutionStep (NOT
        ExecuteInitialize), and it forces ``constrained = true`` (AddEmptyValue
        + SetBool(True)). This test asserts that through the real production
        chain the analytical Bofang temperature is present at the intended
        initialization point and that the water-covered reservoir nodes have
        their TEMPERATURE DOF fixed, while above-water and non-reservoir nodes
        are untouched.
        """
        stages, fixity = self._run_and_read(
            "case/ProjectParameters_production_wrapper.json", record_fixity=True)

        stage = "after_execute_before_solution_loop"
        self.assertIn(stage, stages, "lifecycle stage not recorded")
        self.assertIsNotNone(fixity, "fixity CSV not recorded")

        face = {n: stages[stage][n] for n in self.RESERVOIR_FACE_NODES}
        checked = self._check_analytical_temperature(
            {stage: face}, stage, self.NODE_Y, self.RESERVOIR_FACE_NODES)
        self.assertGreaterEqual(checked, 4)

        # non-reservoir nodes untouched
        tol = 1.0e-9
        for node_id in [17, 18]:
            self.assertAlmostEqual(stages[stage][node_id], 0.0, delta=tol,
                                   msg="non-reservoir node %d got a temperature" % node_id)

        # fixity: water-covered reservoir nodes (y <= 14) are Fixed; the rest are not
        for node_id, is_fixed in fixity[stage].items():
            expected = 1 if (node_id in self.RESERVOIR_FACE_NODES
                             and self.NODE_Y[node_id] <= 14.0) else 0
            self.assertEqual(is_fixed, expected,
                             msg="node %d fixity mismatch" % node_id)


if __name__ == "__main__":
    import unittest
    unittest.main()
