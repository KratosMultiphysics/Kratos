import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
try:
    import torch  # noqa: F401
    import warp  # noqa: F401 - the accelerated neighbor-search backend
    import physicsnemo.nn.functional  # noqa: F401
    have_accelerated = True
except ImportError:
    have_accelerated = False

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


def _RandomCloud(n=30, seed=0):
    rng = numpy.random.default_rng(seed)
    return rng.random((n, 3))


def _EdgeSet(edge_index):
    return {(int(s), int(r)) for s, r in edge_index.T}


class TestParticleBridge(KratosUnittest.TestCase):
    def test_RadiusGraphBruteForce(self):
        positions = _RandomCloud()
        radius = 0.3
        edge_index, edge_features = particle_bridge.BuildParticleGraphFromPositions(
            positions, Kratos.Parameters('{"type": "radius", "radius": %f, "backend": "numpy"}' % radius))

        # reference: every ordered pair within the radius, no self edges
        expected = set()
        for i in range(len(positions)):
            for j in range(len(positions)):
                if i != j and numpy.linalg.norm(positions[i] - positions[j]) <= radius:
                    expected.add((i, j))
        self.assertEqual(_EdgeSet(edge_index), expected)

        # bidirectional and feature contract: relative position + distance
        self.assertEqual(_EdgeSet(edge_index), {(r, s) for s, r in _EdgeSet(edge_index)})
        self.assertEqual(edge_features.shape, (edge_index.shape[1], 4))
        for column in range(edge_index.shape[1]):
            sender, receiver = edge_index[:, column]
            relative = positions[receiver] - positions[sender]
            self.assertTrue(numpy.allclose(edge_features[column, :3], relative))
            self.assertAlmostEqual(
                edge_features[column, 3], numpy.linalg.norm(relative), places=12)

    def test_KnnGraphBruteForce(self):
        positions = _RandomCloud(20, seed=1)
        edge_index, _ = particle_bridge.BuildParticleGraphFromPositions(
            positions, Kratos.Parameters('{"type": "knn", "max_neighbors": 3, "backend": "numpy"}'))
        edges = _EdgeSet(edge_index)
        # symmetrized: every node has >= 3 in-neighbors, no self edges
        for i in range(len(positions)):
            self.assertGreaterEqual(sum(1 for s, r in edges if r == i), 3)
        self.assertTrue(all(s != r for s, r in edges))

    @KratosUnittest.skipUnless(have_accelerated,
                               "Missing required python modules: torch, warp, physicsnemo.")
    def test_AcceleratedMatchesBruteForce(self):
        positions = _RandomCloud(40, seed=2)
        for connectivity in ('{"type": "radius", "radius": 0.35, "backend": "%s"}',
                             '{"type": "knn", "max_neighbors": 4, "backend": "%s"}'):
            accelerated, _ = particle_bridge.BuildParticleGraphFromPositions(
                positions, Kratos.Parameters(connectivity % "auto"))
            brute, _ = particle_bridge.BuildParticleGraphFromPositions(
                positions, Kratos.Parameters(connectivity % "numpy"))
            self.assertEqual(_EdgeSet(accelerated), _EdgeSet(brute))

    def test_ValidationErrors(self):
        positions = _RandomCloud(5)
        with self.assertRaisesRegex(ValueError, "connectivity type"):
            particle_bridge.BuildParticleGraphFromPositions(
                positions, Kratos.Parameters('{"type": "delaunay"}'))
        with self.assertRaisesRegex(ValueError, "radius"):
            particle_bridge.BuildParticleGraphFromPositions(
                positions, Kratos.Parameters('{"type": "radius", "radius": 0.0}'))
        with self.assertRaisesRegex(ValueError, "backend"):
            particle_bridge.BuildParticleGraphFromPositions(
                positions, Kratos.Parameters('{"backend": "octree"}'))

    def test_ModelPartGraphAndKinematicFeatures(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Particles")
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.SetBufferSize(2)
        for i, (x, y, z) in enumerate([(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.9, 0.9, 0.9)]):
            node = model_part.CreateNewNode(i + 1, x, y, z)
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, [float(i + 1), 0.0, 0.0])
            node.SetSolutionStepValue(Kratos.VELOCITY, 1, [10.0 * (i + 1), 0.0, 0.0])

        node_features, edge_index, edge_features, node_ids = particle_bridge.BuildParticleGraph(
            model_part, Kratos.Parameters('{"type": "radius", "radius": 0.2, "backend": "numpy"}'))
        self.assertEqual(list(node_ids), [1, 2, 3])
        self.assertEqual(_EdgeSet(edge_index), {(0, 1), (1, 0)})  # only the close pair
        self.assertEqual(node_features.shape, (3, 0))

        features = particle_bridge.BuildKinematicFeatures(model_part, 2)
        # oldest first: buffer index 1 (10, 20, 30) then index 0 (1, 2, 3)
        expected = numpy.array([
            [10.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [20.0, 0.0, 0.0, 2.0, 0.0, 0.0],
            [30.0, 0.0, 0.0, 3.0, 0.0, 0.0]])
        self.assertTrue(numpy.allclose(features, expected))

        with self.assertRaisesRegex(ValueError, "buffer size"):
            particle_bridge.BuildKinematicFeatures(model_part, 3)


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestParticleTrajectoryDataset(KratosUnittest.TestCase):
    def test_FreeFallWindows(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateParticleTrajectoryDataset

        # free fall: x(t) = x0 - 0.5 g t^2 in z -> constant acceleration -g
        g = 9.81
        dt = 0.1
        steps, n = 8, 4
        rng = numpy.random.default_rng(0)
        x0 = rng.random((n, 3))
        trajectory = numpy.stack([
            x0 - numpy.array([0.0, 0.0, 0.5 * g * (t * dt) ** 2]) for t in range(steps)])

        dataset = CreateParticleTrajectoryDataset(trajectory, history_size=2, delta_time=dt)
        self.assertEqual(len(dataset), steps - 3)  # t in [2, steps-2]
        features, target = dataset[0]
        self.assertEqual(list(features.shape), [n, 6])
        self.assertEqual(list(target.shape), [n, 3])
        # central differences recover the exact constant acceleration
        self.assertTrue(numpy.allclose(target.numpy()[:, 2], -g, atol=1e-9))
        self.assertTrue(numpy.allclose(target.numpy()[:, :2], 0.0, atol=1e-9))
        self.assertEqual(len(dataset.positions), len(dataset))
        self.assertEqual(dataset.target_mean.shape, (3,))

    def test_NormalizationCardEntriesCarryBothHalves(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import (
            CreateParticleTrajectoryDataset, MakeNormalizationCardEntries)
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry

        rng = numpy.random.default_rng(3)
        trajectory = numpy.cumsum(rng.random((8, 4, 3)), axis=0)
        dataset = CreateParticleTrajectoryDataset(trajectory, history_size=2, delta_time=0.1,
                                                  normalize=True)
        entries = MakeNormalizationCardEntries(dataset)
        self.assertEqual(sorted(entries), ["input_normalization", "output_normalization"])
        numpy.testing.assert_allclose(entries["input_normalization"]["mean"], dataset.feature_mean)
        numpy.testing.assert_allclose(entries["input_normalization"]["std"], dataset.feature_std)
        numpy.testing.assert_allclose(entries["output_normalization"]["std"], dataset.target_std)
        # the entries undo exactly what normalize=True did: the first raw
        # window, standardized by the card, is the dataset's first sample
        raw = CreateParticleTrajectoryDataset(trajectory, history_size=2, delta_time=0.1)
        numpy.testing.assert_allclose(
            model_registry.ApplyInputNormalization(raw[0][0].numpy(), entries["input_normalization"]),
            dataset[0][0].numpy(), rtol=1e-10)
        numpy.testing.assert_allclose(
            model_registry.ApplyOutputNormalization(dataset[0][1].numpy(), entries["output_normalization"]),
            raw[0][1].numpy(), rtol=1e-10)
        with self.assertRaisesRegex(ValueError, "feature_mean"):
            MakeNormalizationCardEntries(object())

    def test_Validation(self):
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateParticleTrajectoryDataset

        with self.assertRaisesRegex(ValueError, "history_size"):
            CreateParticleTrajectoryDataset(numpy.zeros((5, 2, 3)), 0, 0.1)
        with self.assertRaisesRegex(ValueError, "delta_time"):
            CreateParticleTrajectoryDataset(numpy.zeros((5, 2, 3)), 2, 0.0)
        with self.assertRaisesRegex(ValueError, "at least"):
            CreateParticleTrajectoryDataset(numpy.zeros((3, 2, 3)), 2, 0.1)
        with self.assertRaisesRegex(ValueError, r"\(T, N, 3\)"):
            CreateParticleTrajectoryDataset(numpy.zeros((5, 2, 2)), 2, 0.1)


if __name__ == '__main__':
    KratosUnittest.main()
