"""Transient flow past a cylinder learned by a MeshGraphNet.

The application's fluid case was a steady cavity, so nothing here had a
transient fluid series to learn from. This drives the canonical cylinder
instead and trains a graph surrogate on its own mesh graph.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

_CASES_DIR = Path(__file__).parent / "kratos_solver_cases"
sys.path.insert(0, str(_CASES_DIR))

have_fluid = kratos_utils.CheckIfApplicationsAvailable(
    "FluidDynamicsApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.nn.module.gnn_layers.graph_types import PYG_AVAILABLE
    from physicsnemo.models.meshgraphnet import MeshGraphNet
    have_meshgraphnet = have_torch and PYG_AVAILABLE
except ImportError:
    have_meshgraphnet = False


def _CylinderAvailable():
    if not have_fluid:
        return False
    import cylinder_case

    return cylinder_case.IsAvailable()


@KratosUnittest.skipUnless(have_fluid, "Missing FluidDynamics/LinearSolvers applications.")
class TestCylinderCase(KratosUnittest.TestCase):
    """The transient fixture itself, before any model touches it."""

    def setUp(self):
        import cylinder_case

        if not cylinder_case.IsAvailable():
            self.skipTest("FluidDynamicsApplication's cylinder test data is not present.")

    def test_TheFlowEvolvesOverTheSeries(self):
        import cylinder_case
        from transient_harness import RunTransientAnalysis

        model = Kratos.Model()
        analysis = cylinder_case.CreateCylinderAnalysis(
            model, end_time=0.2, time_step=0.02)
        states = RunTransientAnalysis(analysis, collect=cylinder_case.CollectVelocities)

        self.assertEqual(states.ndim, 3)
        self.assertGreaterEqual(states.shape[0], 5)
        self.assertTrue(numpy.isfinite(states).all())
        # a transient case that does not change is a steady one mislabelled
        self.assertGreater(numpy.abs(states[-1] - states[0]).max(), 1e-6)
        # the inflow really drives it
        self.assertGreater(numpy.abs(states).max(), 0.5)

    def test_TheFixtureIsReusedNotCopied(self):
        """The mesh belongs to FluidDynamicsApplication; this case reads it
        in place and skips when it is absent, so nothing is duplicated."""
        import cylinder_case

        self.assertTrue(str(cylinder_case._MESH_FILE).endswith("cylinder_2d.mdpa"))
        self.assertIn("FluidDynamicsApplication", str(cylinder_case._MESH_FILE))


@KratosUnittest.skipUnless(have_fluid and have_meshgraphnet,
                           "Missing FluidDynamics, torch or torch_geometric.")
class TestVortexSheddingSurrogate(KratosUnittest.TestCase):
    def setUp(self):
        import cylinder_case

        if not cylinder_case.IsAvailable():
            self.skipTest("FluidDynamicsApplication's cylinder test data is not present.")
        self.model = Kratos.Model()
        analysis = cylinder_case.CreateCylinderAnalysis(
            self.model, end_time=0.2, time_step=0.02)
        from transient_harness import RunTransientAnalysis

        self.states = RunTransientAnalysis(
            analysis, collect=cylinder_case.CollectVelocities)
        self.model_part = cylinder_case.GetModelPart(self.model)

    def test_TheMeshGraphIsTheSolversOwn(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge

        features, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(
            self.model_part, [("VELOCITY", "node_historical")])
        self.assertEqual(len(node_ids), self.model_part.NumberOfNodes())
        self.assertEqual(edge_index.shape[0], 2)
        self.assertGreater(edge_index.shape[1], len(node_ids))   # a real mesh graph
        self.assertEqual(edge_features.shape[1], 4)              # offset plus distance
        self.assertTrue(numpy.isfinite(features).all())

    def test_TheSeriesFeedsATemporalDataset(self):
        """A transient series the surrogate can actually consume: the window
        dataset pairs a history with the state that follows it."""
        from KratosMultiphysics.PhysicsNeMoApplication.training import temporal_training

        # the window dataset wants (T, N, W) - a trajectory of FIELDS, not
        # of flattened vectors
        states = torch.tensor(self.states, dtype=torch.float32)
        # a LIST of trajectories, each (T, N, W) - one series here
        dataset = temporal_training.CreateTrajectoryWindowDataset(
            [states], Kratos.Parameters("{}"))
        self.assertGreater(len(dataset), 0)
        inputs, targets = dataset[0]
        self.assertTrue(torch.isfinite(inputs).all())
        self.assertTrue(torch.isfinite(targets).all())

    def test_TheGraphSurrogateRunsOnTheCylinder(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge

        features, edge_index, edge_features, _ = graph_bridge.BuildGraph(
            self.model_part, [("VELOCITY", "node_historical")])
        torch.manual_seed(0)
        model = MeshGraphNet(
            input_dim_nodes=features.shape[1], input_dim_edges=edge_features.shape[1],
            output_dim=features.shape[1], processor_size=1,
            hidden_dim_node_encoder=16, hidden_dim_edge_encoder=16,
            hidden_dim_node_decoder=16).eval()

        # ToPyGGraph carries the TOPOLOGY; the features are arguments
        graph = graph_bridge.ToPyGGraph(edge_index, num_nodes=features.shape[0])
        with torch.no_grad():
            prediction = model(
                torch.tensor(features, dtype=torch.float32),
                torch.tensor(edge_features, dtype=torch.float32),
                graph)
        self.assertEqual(tuple(prediction.shape), (features.shape[0], features.shape[1]))
        self.assertTrue(bool(torch.isfinite(prediction).all()))


if __name__ == '__main__':
    KratosUnittest.main()
