"""Tests for the Lennard-Jones GNN recipe: the numpy reference integrator
(always runs), the periodic particle graph (always runs), and the
MeshGraphNet force/potential training + ParticleInferenceProcess deployment
(self-skipping without torch_geometric/torch_scatter).

NVIDIA's own examples/molecular_dynamics/lennard_jones trains the generic
MeshGraphNet on (positions -> forces) frames: node feature ones, edge
features the minimum-image relative position plus distance, a radius graph
in a periodic box. No "molecular" architecture exists upstream, and none
is needed - the roadmap had that gate wrong for several rounds. What the
bridge lacked was the periodic neighbour search.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities import lennard_jones_reference as lj

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    import torch_geometric  # noqa: F401
    import torch_scatter  # noqa: F401
    from physicsnemo.models.meshgraphnet import MeshGraphNet
    have_mgn = have_torch
except ImportError:
    have_mgn = False


def _EdgeSet(edge_index):
    return {(int(s), int(r)) for s, r in edge_index.T}


class TestLennardJonesReference(KratosUnittest.TestCase):
    def test_ForceIsZeroAtTheWellMinimumAndAntisymmetric(self):
        box = [10.0]
        r_min = 2.0 ** (1.0 / 6.0)
        forces, potential = lj.ComputeForcesAndPotential(
            [[1.0, 1.0, 1.0], [1.0 + r_min, 1.0, 1.0]], box)
        numpy.testing.assert_allclose(forces, 0.0, atol=1e-12)
        # per-atom energies share the pair energy: each is half of V(r_min) - V(r_c)
        shift = 4.0 * (2.5 ** -12 - 2.5 ** -6)
        numpy.testing.assert_allclose(potential, 0.5 * (-1.0 - shift), rtol=1e-12)
        # closer than the minimum: repulsive, equal and opposite
        forces, _ = lj.ComputeForcesAndPotential([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], box)
        self.assertLess(forces[0, 0], 0.0)
        numpy.testing.assert_allclose(forces[0], -forces[1], rtol=1e-12)

    def test_MinimumImageActsAcrossTheBoundary(self):
        # two atoms 0.2 apart THROUGH the boundary of a box of 5: the
        # direct distance is 4.8, outside the cutoff, so the force is
        # entirely the periodic image's
        box = [5.0]
        forces, _ = lj.ComputeForcesAndPotential([[0.1, 2.0, 2.0], [4.9, 2.0, 2.0]], box)
        self.assertGreater(numpy.abs(forces).max(), 0.0)
        direct, _ = lj.ComputeForcesAndPotential([[0.1, 2.0, 2.0], [0.3, 2.0, 2.0]], box)
        # the image pair (1.0 -> 4.9 is -0.2 away) is the mirror of the direct pair
        numpy.testing.assert_allclose(forces[0], -direct[0], rtol=1e-12)
        numpy.testing.assert_allclose(
            lj.MinimumImage([[4.8, 0.0, 0.0]], box), [[-0.2, 0.0, 0.0]], atol=1e-12)

    def test_EnergyAndMomentumAreConserved(self):
        trajectory = lj.GenerateTrajectory(atoms_per_side=3, steps=60, dt=0.005, seed=1)
        box = trajectory["box_size"]
        energies = [lj.ComputeEnergy(p, v, box)
                    for p, v in zip(trajectory["positions"], trajectory["velocities"])]
        self.assertLess(abs(energies[-1] - energies[0]), 1e-3 * abs(energies[0]))
        momenta = trajectory["velocities"].sum(axis=1)  # mass 1
        numpy.testing.assert_allclose(momenta, 0.0, atol=1e-10)
        # the atoms did move, and the forces are the recorded ones
        self.assertGreater(numpy.abs(trajectory["positions"][-1] - trajectory["positions"][0]).max(), 0.01)
        forces, _ = lj.ComputeForcesAndPotential(trajectory["positions"][10], box)
        numpy.testing.assert_allclose(trajectory["forces"][10], forces, rtol=1e-12)

    def test_CentralDifferencesRecoverTheForcesExactly(self):
        # velocity Verlet IS the Stoermer-Verlet recurrence, so the dataset's
        # acceleration targets are the forces at x_t to round-off
        trajectory = lj.GenerateTrajectory(atoms_per_side=3, steps=6, dt=0.01, seed=2)
        x, dt = trajectory["positions"], trajectory["dt"]
        for t in range(1, 5):
            acceleration = (x[t + 1] - 2.0 * x[t] + x[t - 1]) / dt ** 2
            numpy.testing.assert_allclose(acceleration, trajectory["forces"][t],
                                          rtol=1e-6, atol=1e-8)

    def test_DeterminismAndShapes(self):
        a = lj.GenerateTrajectory(atoms_per_side=2, steps=5, seed=3)
        b = lj.GenerateTrajectory(atoms_per_side=2, steps=5, seed=3)
        numpy.testing.assert_array_equal(a["positions"], b["positions"])
        self.assertEqual(a["positions"].shape, (5, 8, 3))
        self.assertEqual(a["potential"].shape, (5, 8))
        numpy.testing.assert_allclose(a["box_size"], [3.0, 3.0, 3.0])

    def test_Guards(self):
        with self.assertRaisesRegex(ValueError, "cutoff"):
            lj.ComputeForcesAndPotential(numpy.zeros((2, 3)), [4.0], cutoff=2.5)
        with self.assertRaisesRegex(ValueError, "box_size"):
            lj.MinimumImage(numpy.zeros((1, 3)), [1.0, -1.0, 1.0])
        with self.assertRaisesRegex(ValueError, "steps"):
            lj.GenerateTrajectory(steps=1)
        with self.assertRaisesRegex(ValueError, "dt"):
            lj.Step(numpy.zeros((1, 3)), numpy.zeros((1, 3)), 0.0, [5.0])


class TestPeriodicParticleGraph(KratosUnittest.TestCase):
    """"box_size" in the connectivity block: the minimum-image neighbour
    search the molecular-dynamics recipe needs."""

    _STRADDLING = numpy.array([[0.05, 0.5, 0.5], [0.95, 0.5, 0.5], [0.5, 0.5, 0.5]])

    def test_APairStraddlingTheBoundaryIsLinkedWithTheShortVector(self):
        for backend in ("numpy", "auto"):
            with self.subTest(backend=backend):
                edge_index, edge_features = particle_bridge.BuildParticleGraphFromPositions(
                    self._STRADDLING, Kratos.Parameters(
                        '{"type": "radius", "radius": 0.2, "box_size": [1.0], "backend": "%s"}' % backend))
                self.assertEqual(_EdgeSet(edge_index), {(0, 1), (1, 0)})
                for column in range(edge_index.shape[1]):
                    sender, receiver = edge_index[:, column]
                    expected = -0.1 if (sender, receiver) == (0, 1) else 0.1
                    numpy.testing.assert_allclose(edge_features[column, :3], [expected, 0.0, 0.0],
                                                  atol=1e-12)
                    self.assertAlmostEqual(edge_features[column, 3], 0.1, places=12)

    def test_WithoutABoxTheSamePairIsNotLinked(self):
        edge_index, _ = particle_bridge.BuildParticleGraphFromPositions(
            self._STRADDLING, Kratos.Parameters('{"type": "radius", "radius": 0.2, "backend": "numpy"}'))
        self.assertEqual(edge_index.shape, (2, 0))

    def test_UnwrappedPositionsAreTakenModuloTheBox(self):
        shifted = self._STRADDLING + numpy.array([[3.0, -2.0, 0.0]])  # outside [0, 1)
        edge_index, edge_features = particle_bridge.BuildParticleGraphFromPositions(
            shifted, Kratos.Parameters('{"type": "radius", "radius": 0.2, "box_size": [1.0], "backend": "auto"}'))
        self.assertEqual(_EdgeSet(edge_index), {(0, 1), (1, 0)})
        self.assertAlmostEqual(edge_features[0, 3], 0.1, places=12)

    def test_TreeMatchesBruteForceOnARandomPeriodicCloud(self):
        rng = numpy.random.default_rng(0)
        positions = rng.random((60, 3)) * numpy.array([2.0, 3.0, 1.5])
        for connectivity in ('{"type": "radius", "radius": 0.45, "box_size": [2.0, 3.0, 1.5], "backend": "%s"}',
                             '{"type": "knn", "max_neighbors": 4, "box_size": [2.0, 3.0, 1.5], "backend": "%s"}'):
            with self.subTest(connectivity=connectivity[:16]):
                tree, tree_features = particle_bridge.BuildParticleGraphFromPositions(
                    positions, Kratos.Parameters(connectivity % "auto"))
                brute, brute_features = particle_bridge.BuildParticleGraphFromPositions(
                    positions, Kratos.Parameters(connectivity % "numpy"))
                self.assertEqual(_EdgeSet(tree), _EdgeSet(brute))
                numpy.testing.assert_allclose(tree_features, brute_features, atol=1e-12)
                self.assertGreater(tree.shape[1], 0)
        # and the brute-force edge vectors really are minimum images
        edge_index, edge_features = particle_bridge.BuildParticleGraphFromPositions(
            positions, Kratos.Parameters(
                '{"type": "radius", "radius": 0.45, "box_size": [2.0, 3.0, 1.5], "backend": "numpy"}'))
        self.assertLessEqual(numpy.abs(edge_features[:, :3]).max(), 0.45 + 1e-12)

    def test_Guards(self):
        with self.assertRaisesRegex(ValueError, "box_size"):
            particle_bridge.BuildParticleGraphFromPositions(
                self._STRADDLING, Kratos.Parameters('{"box_size": [1.0, 1.0]}'))
        with self.assertRaisesRegex(ValueError, "half the box"):
            particle_bridge.BuildParticleGraphFromPositions(
                self._STRADDLING, Kratos.Parameters('{"type": "radius", "radius": 0.6, "box_size": [1.0]}'))


@KratosUnittest.skipUnless(have_mgn,
                           "Missing required python modules: physicsnemo with torch_geometric/torch_scatter.")
class TestLennardJonesMeshGraphNet(KratosUnittest.TestCase):
    """The recipe end to end: MeshGraphNet trained on periodic radius graphs
    of reference frames to predict per-atom acceleration (force, mass 1)
    and per-atom potential energy, saved with a two-key card, deployed
    through ParticleInferenceProcess with a periodic box."""

    _CONNECTIVITY = '{"type": "radius", "radius": 2.5, "box_size": [6.0], "backend": "auto"}'

    def setUp(self):
        self.checkpoint = Path("test_lennard_jones_mgn.mdlus")
        self.trajectory = lj.GenerateTrajectory(atoms_per_side=4, steps=14, dt=0.005, seed=0)
        self.dt = self.trajectory["dt"]

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint) + ".card.json")

    @staticmethod
    def _Model(output_dim, seed=0):
        torch.manual_seed(seed)
        return MeshGraphNet(input_dim_nodes=6, input_dim_edges=4, output_dim=output_dim,
                            processor_size=2, hidden_dim_processor=16,
                            hidden_dim_node_encoder=16, hidden_dim_edge_encoder=16,
                            hidden_dim_node_decoder=16).double()

    def _Samples(self, dataset):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
        samples = []
        for i in range(len(dataset)):
            features, target = dataset[i]
            edge_index, edge_features = particle_bridge.BuildParticleGraphFromPositions(
                dataset.positions[i], Kratos.Parameters(self._CONNECTIVITY))
            samples.append((features, torch.from_numpy(edge_features),
                            graph_bridge.ToPyGGraph(edge_index, features.shape[0]), target))
        return samples

    @staticmethod
    def _Train(model, samples, epochs, lr=2e-3):
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        history = []
        for _ in range(epochs):
            total = 0.0
            for features, edge_features, graph, target in samples:
                optimizer.zero_grad()
                loss = torch.nn.functional.mse_loss(model(features, edge_features, graph), target)
                loss.backward()
                optimizer.step()
                total += loss.item()
            history.append(total / len(samples))
        return history

    def test_TrainForcesSaveWithCardAndDeployPeriodically(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import (
            CreateParticleTrajectoryDataset, MakeNormalizationCardEntries)
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import particle_inference_process

        dataset = CreateParticleTrajectoryDataset(
            self.trajectory["positions"], history_size=2, delta_time=self.dt, normalize=True)
        # the targets ARE the reference forces (velocity Verlet is exact here)
        raw = CreateParticleTrajectoryDataset(
            self.trajectory["positions"], history_size=2, delta_time=self.dt)
        numpy.testing.assert_allclose(raw[0][1].numpy(), self.trajectory["forces"][2],
                                      rtol=1e-6, atol=1e-8)

        samples = self._Samples(dataset)
        self.assertTrue(all(s[2].edge_index.shape[1] > 0 for s in samples))
        model = self._Model(output_dim=3)
        history = self._Train(model, samples, epochs=30)
        self.assertLess(history[-1], 0.8 * history[0])

        card = {"input_fields": [{"variable_name": "VELOCITY", "data_location": "node_historical"}],
                "output_fields": [{"variable_name": "ACCELERATION", "data_location": "node_historical"}],
                "history_size": 2, "connectivity": {"radius": 2.5, "box_size": [6.0]}}
        card.update(MakeNormalizationCardEntries(dataset))
        training_utils.SaveTrainedModel(model, self.checkpoint, card=card)

        # deployment on a Kratos particle cloud at frame t, with the true
        # window (v_{t-1}, v_t) fed over two steps
        t = 5
        x, v = self.trajectory["positions"], self.trajectory["velocities"]
        model_kratos = Kratos.Model()
        model_part = model_kratos.CreateModelPart("Atoms")
        for variable in (Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.DISPLACEMENT):
            model_part.AddNodalSolutionStepVariable(variable)
        model_part.SetBufferSize(2)
        for i, xyz in enumerate(x[t]):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.ProcessInfo[Kratos.DELTA_TIME] = self.dt
        process = particle_inference_process.Factory(Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Atoms",
                "model_settings"  : { "checkpoint_file" : "%s", "checkpoint_type" : "physicsnemo",
                                      "device" : "cpu" },
                "model_interface" : "meshgraphnet",
                "connectivity"    : %s,
                "history_size"    : 2
            }
        }""" % (self.checkpoint, self._CONNECTIVITY)), model_kratos)
        finite_velocity = [(x[t - 1] - x[t - 2]) / self.dt, (x[t] - x[t - 1]) / self.dt]
        for step, velocity in enumerate(finite_velocity, start=1):
            for i, node in enumerate(model_part.Nodes):
                node.SetSolutionStepValue(Kratos.VELOCITY, [float(c) for c in velocity[i]])
            model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()

        predicted = numpy.array([node.GetSolutionStepValue(Kratos.ACCELERATION)
                                 for node in model_part.Nodes])
        reference = self.trajectory["forces"][t]
        self.assertTrue(numpy.isfinite(predicted).all())
        error = numpy.sqrt(numpy.mean((predicted - reference) ** 2))
        # better than predicting the mean force everywhere - a model that
        # ignored its (de-normalized) inputs could not do this
        baseline = numpy.sqrt(numpy.mean((reference - reference.mean(axis=0)) ** 2))
        self.assertLess(error, baseline)
        # and the nodes advanced (unwrapped; the periodic graph does not care)
        self.assertGreater(numpy.abs(numpy.array([node.X for node in model_part.Nodes]) - x[t][:, 0]).max(), 0.0)

    def test_PotentialEnergyHead(self):
        # the second half of "force and potential prediction": per-atom
        # potential energies as a graph-node regression on the same graphs
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import (
            CreateParticleTrajectoryDataset)
        dataset = CreateParticleTrajectoryDataset(
            self.trajectory["positions"], history_size=2, delta_time=self.dt, normalize=True)
        potential = self.trajectory["potential"][2:2 + len(dataset)]
        mean, std = potential.mean(), potential.std()
        samples = [(f, e, g, torch.from_numpy((potential[i] - mean) / std)[:, None])
                   for i, (f, e, g, _) in enumerate(self._Samples(dataset))]
        model = self._Model(output_dim=1, seed=1)
        history = self._Train(model, samples, epochs=15)
        self.assertLess(history[-1], 0.8 * history[0])


if __name__ == '__main__':
    KratosUnittest.main()
