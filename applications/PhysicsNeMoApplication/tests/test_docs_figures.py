"""The documentation figure generator must keep producing its files.

Runs the numpy-and-matplotlib subset of examples/figures/make_docs_figures.py
into a temporary directory (the Kratos-driven figures are exercised by the
--list contract and by the bridge tests they call into). Skips cleanly where
matplotlib is not installed; never needs torch.
"""

import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest

_SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "figures" / "make_docs_figures.py"
_HAVE_MATPLOTLIB = importlib.util.find_spec("matplotlib") is not None


@KratosUnittest.skipUnless(_HAVE_MATPLOTLIB, "matplotlib is not available")
class TestDocsFigures(KratosUnittest.TestCase):

    def test_ListNamesEveryFigure(self):
        completed = subprocess.run([sys.executable, str(_SCRIPT), "--list"], capture_output=True, text=True, timeout=300)
        self.assertEqual(completed.returncode, 0, completed.stderr[-2000:])
        for name in ("tessellation_modes", "higher_order_modes", "mesh_graph", "proximity_graph",
                     "grid_sampling", "calibration_views", "halo_partition", "benchmark_costs"):
            self.assertIn(name, completed.stdout)

    def test_TorchFreeFiguresAreGenerated(self):
        with tempfile.TemporaryDirectory() as out_dir:
            completed = subprocess.run(
                [sys.executable, str(_SCRIPT), "--out", out_dir, "--only", "calibration_views,halo_partition,proximity_graph,higher_order_modes"],
                capture_output=True, text=True, timeout=900, env={**os.environ, "MPLBACKEND": "Agg"})
            self.assertEqual(completed.returncode, 0, completed.stdout[-2000:] + completed.stderr[-2000:])
            for relative in ("Uncertainty/images/calibration_views.png", "Distributed/images/halo_partition.png",
                             "Particle_Methods/images/proximity_graph.png", "Mesh_Bridge/images/higher_order_modes.png"):
                path = Path(out_dir) / relative
                self.assertTrue(path.is_file(), f"{relative} was not produced")
                self.assertGreater(path.stat().st_size, 10_000, f"{relative} is suspiciously small")


if __name__ == "__main__":
    KratosUnittest.main()
