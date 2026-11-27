"""Tests for the GraphCast shallow-water recipe: the numpy reference
integrator (always runs) and the GraphCastNet training/deployment path
(self-skipping - GraphCastNet needs torch_geometric PLUS torch_sparse, or
dgl, at CONSTRUCTION time)."""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.utilities import shallow_water_reference

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

# GraphCastNet's graph backend needs torch_geometric + torch_sparse (or
# dgl); the ImportError fires at CONSTRUCTION, so probe with a tiny build.
have_graphcast = False
if have_torch:
    try:
        from physicsnemo.models.graphcast import GraphCastNet
        GraphCastNet(mesh_level=0, input_res=(4, 8),
                     input_dim_grid_nodes=1, output_dim_grid_nodes=1,
                     processor_layers=1, hidden_dim=4, hidden_layers=1)
        have_graphcast = True
    except Exception:
        have_graphcast = False


class TestShallowWaterReference(KratosUnittest.TestCase):
    def test_TrajectoryShapesAndDeterminism(self):
        trajectory = shallow_water_reference.GenerateTrajectory(
            shape=(8, 16), steps=10, seed=3)
        self.assertEqual(trajectory.shape, (10, 3, 8, 16))
        again = shallow_water_reference.GenerateTrajectory(
            shape=(8, 16), steps=10, seed=3)
        numpy.testing.assert_array_equal(trajectory, again)

        pairs = shallow_water_reference.MakeStepPairs(trajectory)
        self.assertEqual(len(pairs), 9)
        numpy.testing.assert_array_equal(pairs[0][1], pairs[1][0])

    def test_FlatStateIsAFixedPoint(self):
        flat = numpy.zeros((3, 6, 6))
        stepped = shallow_water_reference.Step(
            flat, shallow_water_reference.StableTimeStep())
        numpy.testing.assert_array_equal(stepped, flat)

    def test_EnergyStaysBounded(self):
        trajectory = shallow_water_reference.GenerateTrajectory(
            shape=(12, 12), steps=60, seed=0)
        energies = [shallow_water_reference.ComputeEnergy(state) for state in trajectory]
        self.assertGreater(energies[0], 0.0)
        # RK2 under the CFL guard: the quadratic invariant drifts only weakly
        self.assertLess(max(energies), 1.05 * energies[0])
        self.assertGreater(min(energies), 0.5 * energies[0])

    def test_MassIsConserved(self):
        # periodic domain: the height integral is exactly conserved by the
        # centered divergence (roll sums telescope)
        trajectory = shallow_water_reference.GenerateTrajectory(
            shape=(8, 8), steps=30, seed=1)
        masses = [float(state[0].sum()) for state in trajectory]
        numpy.testing.assert_allclose(masses, masses[0], atol=1e-10)

    def test_CflGuard(self):
        with self.assertRaisesRegex(ValueError, "CFL"):
            shallow_water_reference.GenerateTrajectory(shape=(8, 8), steps=5, dt=10.0)
        with self.assertRaisesRegex(ValueError, "steps"):
            shallow_water_reference.GenerateTrajectory(steps=1)
        with self.assertRaisesRegex(ValueError, "initial_state"):
            shallow_water_reference.GenerateTrajectory(
                shape=(8, 8), steps=5, initial_state=numpy.zeros((3, 4, 4)))


@KratosUnittest.skipUnless(have_graphcast,
                           "GraphCastNet unavailable (needs torch + physicsnemo + "
                           "torch_geometric with torch_sparse, or dgl).")
class TestGraphCastShallowWater(KratosUnittest.TestCase):
    """The full recipe: tiny GraphCastNet trained on step pairs (B=1!),
    saved as .mdlus, deployed through GridInferenceProcess's "grid"
    interface with the squeeze idiom."""

    def setUp(self):
        self.checkpoint = Path("test_graphcast_sw.mdlus")

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint) + ".card.json")

    def _TinyGraphCast(self):
        from physicsnemo.models.graphcast import GraphCastNet

        torch.manual_seed(0)
        return GraphCastNet(
            mesh_level=1, input_res=(8, 16),
            input_dim_grid_nodes=3, output_dim_grid_nodes=3,
            processor_layers=3, hidden_dim=16, hidden_layers=1)

    def test_TrainSaveDeploy(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import grid_inference_process
        from test_grid_bridge import CreateStructuredTetModelPart

        model = self._TinyGraphCast()
        with torch.no_grad():
            out = model(torch.zeros(1, 3, 8, 16))
        self.assertEqual(tuple(out.shape), (1, 3, 8, 16))

        trajectory = shallow_water_reference.GenerateTrajectory(
            shape=(8, 16), steps=8, seed=0)
        pairs = shallow_water_reference.MakeStepPairs(trajectory)
        inputs = torch.tensor(numpy.stack([p[0] for p in pairs]), dtype=torch.float32)
        targets = torch.tensor(numpy.stack([p[1] for p in pairs]), dtype=torch.float32)
        dataset = torch.utils.data.TensorDataset(inputs, targets)

        # GraphCastNet's forward is batch-size-1 only: batch_size MUST be 1
        history = training_utils.TrainModel(model, dataset, Kratos.Parameters("""{
            "epochs"        : 3,
            "batch_size"    : 1,
            "learning_rate" : 1e-3,
            "device"        : "cpu",
            "seed"          : 0
        }"""))
        self.assertLess(history[-1], history[0])
        training_utils.SaveTrainedModel(model, self.checkpoint)

        # deployment: the squeezed (1, 3, 8, 16) batch of the "grid"
        # interface matches GraphCastNet's forward exactly
        kratos_model = Kratos.Model()
        CreateStructuredTetModelPart(
            kratos_model, "Main", divisions=2,
            historical_variables=(Kratos.VELOCITY, Kratos.MESH_VELOCITY))
        model_part = kratos_model["Main"]
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, [node.X, node.Y, 0.1])

        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "%s",
                    "checkpoint_type" : "physicsnemo",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "VELOCITY",      "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "MESH_VELOCITY", "data_location" : "node_historical" } ],
                "grid_shape"      : [8, 16, 2],
                "squeeze_axis"    : 2
            }
        }""")
        process = grid_inference_process.Factory(settings, kratos_model)
        model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        values = numpy.array([
            list(node.GetSolutionStepValue(Kratos.MESH_VELOCITY))
            for node in model_part.Nodes])
        self.assertTrue(numpy.isfinite(values).all())
        self.assertGreater(numpy.abs(values).max(), 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
