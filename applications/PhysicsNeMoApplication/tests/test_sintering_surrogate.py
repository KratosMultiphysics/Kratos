"""The VFGN sintering bridge driven by a REAL thermo-mechanical solve.

The bridge's own tests use synthetic shrinking clouds; this exercises it on
trajectories produced by an actual coupled solver, where the shrinkage comes
out of the thermal strain rather than being prescribed.
"""

import sys
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset
from KratosMultiphysics.PhysicsNeMoApplication.bridges import vfgn_bridge
sys.path.insert(0, str(Path(__file__).parent / "kratos_solver_cases"))

have_coupled = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "StructuralMechanicsApplication",
    "ConstitutiveLawsApplication", "LinearSolversApplication")

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import torch_scatter  # noqa: F401
    from physicsnemo.models.vfgn import VFGNLearnedSimulator  # noqa: F401
    have_vfgn = have_torch
except ImportError:
    have_vfgn = False

_MISSING_APPS = ("Missing required applications: ConvectionDiffusion, "
                 "StructuralMechanics, ConstitutiveLaws, LinearSolvers.")
_MISSING_VFGN = "Missing required python modules: torch, torch_scatter, physicsnemo."


def _RunSintering(divisions=6, end_time=0.5):
    import thermomechanical_case
    import transient_harness

    model = Kratos.Model()
    analysis = thermomechanical_case.CreateSinteringAnalysis(
        model, divisions=divisions, end_time=end_time)
    positions = transient_harness.RunTransientAnalysis(
        analysis, collect=thermomechanical_case.CollectPositions)
    return model, positions


def _Span(state):
    return float(state[:, 0].max() - state[:, 0].min())


@KratosUnittest.skipUnless(have_coupled, _MISSING_APPS)
class TestSinteringCase(KratosUnittest.TestCase):
    """The physics has to be real, not merely free of errors."""

    def test_CoolingContractsTheBody(self):
        import thermomechanical_case

        model, positions = _RunSintering()
        self.assertEqual(positions.ndim, 3)
        self.assertEqual(positions.shape[2], 3)
        self.assertGreater(positions.shape[0], 5)

        # the body contracts, and keeps contracting
        spans = [_Span(state) for state in positions]
        self.assertLess(spans[-1], spans[0])
        self.assertGreater((spans[0] - spans[-1]) / spans[0], 0.05)
        for earlier, later in zip(spans, spans[1:]):
            self.assertLessEqual(later, earlier + 1e-9)

        # ... driven by the temperature actually dropping
        model_part = model["Structure"]
        temperatures = thermomechanical_case.CollectTemperatures(model_part)
        self.assertLess(max(temperatures), 1000.0)

    def test_ShrinkageTracksTheCoolingRate(self):
        """A faster ramp must contract more by the same time - the coupling
        is doing the work, not an arbitrary prescribed motion."""
        import thermomechanical_case
        import transient_harness

        contractions = []
        for cooling_rate in (400.0, 1600.0):
            model = Kratos.Model()
            analysis = thermomechanical_case.CreateSinteringAnalysis(
                model, divisions=5, cooling_rate=cooling_rate, end_time=0.3)
            states = transient_harness.RunTransientAnalysis(
                analysis, collect=thermomechanical_case.CollectPositions)
            contractions.append(_Span(states[0]) - _Span(states[-1]))
        self.assertGreater(contractions[1], contractions[0])

    def test_DeformedPositionsDifferFromTheReference(self):
        """move_mesh_flag is on, so node.X is deformed while X0 is not."""
        model, positions = _RunSintering(divisions=5, end_time=0.3)
        model_part = model["Structure"]
        reference = numpy.array([[n.X0, n.Y0, n.Z0] for n in model_part.Nodes])
        current = numpy.array([[n.X, n.Y, n.Z] for n in model_part.Nodes])
        self.assertGreater(float(numpy.abs(current - reference).max()), 1e-6)
        numpy.testing.assert_allclose(current, positions[-1], atol=1e-12)


@KratosUnittest.skipUnless(have_coupled and have_vfgn,
                           _MISSING_APPS + " " + _MISSING_VFGN)
class TestVfgnOnRealSinteringTrajectories(KratosUnittest.TestCase):

    def setUp(self):
        self.model, self.states = _RunSintering(divisions=5, end_time=0.5)
        # VFGN wants node-major (N, T, 3); the harness collects (T, N, 3)
        self.positions = numpy.transpose(self.states, (1, 0, 2))
        # the radius must resolve the element size (~1/divisions); the
        # bridge's 0.015 default would give an empty graph on this mesh
        self.connectivity = Kratos.Parameters(
            '{"type": "radius", "radius": 0.35, "backend": "numpy"}')

    def _Edges(self, state):
        edge_index, _ = particle_bridge.BuildParticleGraphFromPositions(
            state, self.connectivity)
        return edge_index

    def _Stats(self):
        dataset = torch_dataset.CreateParticleTrajectoryDataset(
            self.states, history_size=2, delta_time=0.05)
        return vfgn_bridge.StatsFromTrajectoryDataset(dataset)

    def test_StatsComeFromTheRealTrajectory(self):
        stats = self._Stats()
        self.assertEqual(tuple(numpy.asarray(stats["velocity"].mean).shape), (3,))
        # a contracting body has genuinely non-zero velocity statistics
        self.assertGreater(float(numpy.abs(numpy.asarray(stats["velocity"].mean)).max()), 0.0)

    def test_TrainingStepOnRealTrajectoriesReducesTheLoss(self):
        torch.manual_seed(0)
        model = vfgn_bridge.CreateVfgnSimulator(Kratos.Parameters("""{
            "predict_length"     : 1,
            "num_seq"            : 5,
            "connectivity_param" : 0.35,
            "boundaries"         : [[-0.5, 1.5], [-0.5, 1.5], [-0.5, 1.5]]
        }"""), self._Stats())

        sequence = self.positions[:, :5, :]
        next_positions = self.positions[:, 5:6, :]
        edges = self._Edges(self.states[4])

        first = vfgn_bridge.ComputeVfgnLoss(
            model, sequence, next_positions, edges, predict_length=1)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)  # after the lazy build
        history = [float(first)]
        for _ in range(12):
            optimizer.zero_grad()
            loss = vfgn_bridge.ComputeVfgnLoss(
                model, sequence, next_positions, edges, predict_length=1)
            loss.backward()
            optimizer.step()
            history.append(float(loss))
        self.assertLess(min(history[1:]), history[0])

    def test_RolloutWithTheFurnaceScheduleAsContext(self):
        """The sintering-specific input the synthetic tests cannot exercise:
        the temperature schedule enters as VFGN's global context."""
        import thermomechanical_case

        torch.manual_seed(0)
        model = vfgn_bridge.CreateVfgnSimulator(Kratos.Parameters("""{
            "predict_length"     : 1,
            "num_seq"            : 5,
            "connectivity_param" : 0.35,
            "boundaries"         : [[-0.5, 1.5], [-0.5, 1.5], [-0.5, 1.5]]
        }"""), self._Stats())

        temperatures = thermomechanical_case.CollectTemperatures(self.model["Structure"])
        furnace = numpy.array([[float(numpy.mean(temperatures)) / 1000.0]])

        predicted = vfgn_bridge.RunVfgnRollout(
            model, self.positions[:, :5, :], self._Edges(self.states[4]),
            predict_length=1, global_context=furnace)
        self.assertEqual(predicted.shape, (self.positions.shape[0], 1, 3))
        self.assertTrue(numpy.isfinite(predicted).all())

    def test_GraphIsNonEmptyAtTheChosenRadius(self):
        """Guards the trap that motivated overriding connectivity_param."""
        edges = self._Edges(self.states[0])
        self.assertGreater(edges.shape[1], 0)
        default_edges, _ = particle_bridge.BuildParticleGraphFromPositions(
            self.states[0], Kratos.Parameters(
                '{"type": "radius", "radius": 0.015, "backend": "numpy"}'))
        self.assertEqual(default_edges.shape[1], 0)   # the default would be empty


if __name__ == '__main__':
    KratosUnittest.main()
