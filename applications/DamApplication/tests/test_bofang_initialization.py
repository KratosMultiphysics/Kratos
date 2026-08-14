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

    def _run_and_read(self):
        with KratosUnittest.WorkFolderScope("bofang_initialization_regression", __file__):
            out_dir = tempfile.mkdtemp(prefix="bofang_regression_")
            RunInstrumentedAnalysis(
                os.path.join("case", "ProjectParameters.json"),
                out_dir,
                "test",
                self.MONITORED_NODES,
                legacy_process_initialize=False,
            )
            lifecycle_path = os.path.join(out_dir, "lifecycle_temperatures.csv")
            self.assertTrue(os.path.exists(lifecycle_path), "lifecycle CSV was not written")
            rows = {}
            with open(lifecycle_path) as f:
                header = f.readline().strip().split(",")
                for line in f:
                    parts = line.strip().split(",")
                    record = dict(zip(header, parts))
                    stage = record["lifecycle_stage"]
                    node_id = int(record["node_id"])
                    rows.setdefault(stage, {})[node_id] = float(record["temperature"])
            return rows

    def test_bofang_initial_temperature_matches_analytical(self):
        stages = self._run_and_read()

        analytical = BofangAnalytical()

        # stage reached right after the whole Initialize() (post-solver-init
        # ExecuteBeforeSolutionLoop), i.e. the state at the first solution step
        stage = "after_execute_before_solution_loop"
        self.assertIn(stage, stages, "lifecycle stage not recorded")

        tol = 1.0e-9
        checked = 0
        for node_id, temp in stages[stage].items():
            if node_id not in self.RESERVOIR_FACE_NODES:
                # not part of the reservoir model part: must stay untouched
                self.assertAlmostEqual(temp, 0.0, delta=tol,
                                       msg="non-reservoir node %d got a temperature" % node_id)
                continue
            expected = analytical.temperature({
                "id": node_id, "x": 0.0, "y": self.NODE_Y[node_id]
            })
            if expected is None:
                # above the water level: the process must not assign anything
                self.assertAlmostEqual(temp, 0.0, delta=tol,
                                       msg="above-water node %d got a temperature" % node_id)
            else:
                self.assertAlmostEqual(temp, expected, delta=tol,
                                       msg="node %d does not match analytical Bofang" % node_id)
                checked += 1

        # several reservoir nodes must genuinely carry the depth-dependent profile
        self.assertGreaterEqual(checked, 4)

    def test_bofang_temperature_available_before_solver_initialize(self):
        stages = self._run_and_read()
        pre = stages.get("before_solver_initialize", {})
        post = stages.get("after_execute_before_solution_loop", {})
        # the water-covered face node (7: y=12 below water level 14) must carry
        # the Bofang value already before the solver is initialized
        self.assertIn(7, pre)
        self.assertAlmostEqual(pre[7], post[7], delta=1.0e-12,
                               msg="Bofang temperature changed during solver initialize")


if __name__ == "__main__":
    import unittest
    unittest.main()
