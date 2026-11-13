"""Every CuPy-converted path must agree with its numpy reference.

numpy is the reference implementation at each of these sites; the GPU path
is an optimization and nothing more, so the contract is agreement to
tolerance (not bit-for-bit - CuPy reorders floating-point reductions, which
is exactly why "auto" never selects it on its own).

The size thresholds are patched to zero here so the device path actually
runs on test-sized inputs; in production those thresholds are what keep a
small mesh on the CPU.
"""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import rom_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils

try:
    import cupy
    have_cupy = array_backend_utils.IsCuPyAvailable()
except ImportError:
    have_cupy = False

_MISSING = "Missing required python module: cupy (with a CUDA device)."


class _ForcedDeviceMixin:
    """Drops both size thresholds so "cupy" really means cupy here.

    Deliberately not a TestCase subclass: it carries no tests, and the
    registration guard rightly refuses TestCase classes that never reach a
    suite.
    """

    def setUp(self):
        super().setUp()
        self._threshold = array_backend_utils.DEFAULT_SIZE_THRESHOLD
        self._rom_threshold = rom_bridge._GPU_BASIS_THRESHOLD
        array_backend_utils.DEFAULT_SIZE_THRESHOLD = 0
        rom_bridge._GPU_BASIS_THRESHOLD = 0

    def tearDown(self):
        array_backend_utils.DEFAULT_SIZE_THRESHOLD = self._threshold
        rom_bridge._GPU_BASIS_THRESHOLD = self._rom_threshold
        super().tearDown()


@KratosUnittest.skipUnless(have_cupy, _MISSING)
class TestCuPyBackendParity(_ForcedDeviceMixin, KratosUnittest.TestCase):

    def test_ParticleRadiusGraphMatchesNumpy(self):
        rng = numpy.random.default_rng(0)
        positions = rng.random((120, 3))
        settings = Kratos.Parameters(
            """{"type" : "radius", "radius" : 0.25, "backend" : "numpy"}""")
        reference_index, reference_features = particle_bridge.BuildParticleGraphFromPositions(
            positions, settings)

        settings = Kratos.Parameters(
            """{"type" : "radius", "radius" : 0.25, "backend" : "cupy"}""")
        index, features = particle_bridge.BuildParticleGraphFromPositions(positions, settings)

        self.assertTrue(numpy.array_equal(reference_index, index))
        self.assertEqual(reference_features.shape, features.shape)
        self.assertTrue(numpy.allclose(reference_features, features, rtol=0.0, atol=1e-12))

    def test_ParticleKnnGraphMatchesNumpy(self):
        rng = numpy.random.default_rng(1)
        positions = rng.random((80, 3))
        reference, _ = particle_bridge.BuildParticleGraphFromPositions(
            positions, Kratos.Parameters(
                """{"type" : "knn", "max_neighbors" : 6, "backend" : "numpy"}"""))
        index, _ = particle_bridge.BuildParticleGraphFromPositions(
            positions, Kratos.Parameters(
                """{"type" : "knn", "max_neighbors" : 6, "backend" : "cupy"}"""))
        self.assertTrue(numpy.array_equal(reference, index))

    def test_RomProjectionRoundTripMatchesNumpy(self):
        rng = numpy.random.default_rng(2)
        phi, _ = numpy.linalg.qr(rng.random((200, 8)))
        basis = rom_bridge.RomBasis(phi=phi, node_ids=numpy.arange(1, 101, dtype=numpy.int64),
                                    nodal_unknowns=("DISPLACEMENT_X", "DISPLACEMENT_Y"),
                                    singular_values=None)
        u = rng.random(200)

        reference_q = rom_bridge.ProjectToReducedSpace(basis, u, backend="numpy")
        device_q = rom_bridge.ProjectToReducedSpace(basis, u, backend="cupy")
        self.assertVectorAlmostEqual(Kratos.Vector(reference_q.tolist()),
                                     Kratos.Vector(device_q.tolist()), 10)

        reference_u = rom_bridge.ReconstructFromReducedSpace(basis, reference_q, backend="numpy")
        device_u = rom_bridge.ReconstructFromReducedSpace(basis, device_q, backend="cupy")
        self.assertVectorAlmostEqual(Kratos.Vector(reference_u.tolist()),
                                     Kratos.Vector(device_u.tolist()), 10)

    def test_TheResidentBasisIsUploadedOnceAndStaysCorrect(self):
        """The device copy is cached on the (frozen) dataclass; repeated
        calls must keep answering with the same values."""
        rng = numpy.random.default_rng(3)
        phi = rng.random((60, 4))
        basis = rom_bridge.RomBasis(phi=phi, node_ids=numpy.arange(1, 31, dtype=numpy.int64),
                                    nodal_unknowns=("VELOCITY_X", "VELOCITY_Y"),
                                    singular_values=None)
        q = rng.random(4)
        reference = rom_bridge.ReconstructFromReducedSpace(basis, q, backend="numpy")
        self.assertIsNone(basis._device_phi)
        for _ in range(3):
            device = rom_bridge.ReconstructFromReducedSpace(basis, q, backend="cupy")
            self.assertVectorAlmostEqual(Kratos.Vector(reference.tolist()),
                                         Kratos.Vector(device.tolist()), 10)
        self.assertIsNotNone(basis._device_phi)

    def test_TrilinearInterpolationMatchesNumpy(self):
        rng = numpy.random.default_rng(4)
        grid = rng.random((3, 9, 8, 7))
        bounding_box = (numpy.zeros(3), numpy.ones(3))
        points = rng.random((250, 3))
        reference = grid_bridge.InterpolateGridAtPoints(
            grid, bounding_box, points, backend="numpy")
        device = grid_bridge.InterpolateGridAtPoints(
            grid, bounding_box, points, backend="cupy")
        self.assertEqual(reference.shape, device.shape)
        self.assertTrue(numpy.allclose(reference, device, rtol=0.0, atol=1e-12))

    def test_TrilinearInterpolationMatchesOutsideTheBox(self):
        """Points are clamped; the clamp must clamp identically on both."""
        rng = numpy.random.default_rng(5)
        grid = rng.random((2, 6, 6, 6))
        bounding_box = (numpy.zeros(3), numpy.ones(3))
        points = rng.random((100, 3)) * 3.0 - 1.0  # well outside
        reference = grid_bridge.InterpolateGridAtPoints(
            grid, bounding_box, points, backend="numpy")
        device = grid_bridge.InterpolateGridAtPoints(
            grid, bounding_box, points, backend="cupy")
        self.assertTrue(numpy.allclose(reference, device, rtol=0.0, atol=1e-12))

    def test_EdgeFeaturesMatchNumpy(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Main")
        properties = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                 (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 1.0)]):
            model_part.CreateNewNode(i + 1, *xyz)
        model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)
        model_part.CreateNewElement("Element3D4N", 2, [2, 5, 3, 4], properties)

        _, edge_index, reference, _ = graph_bridge.BuildGraph(model_part)
        device = graph_bridge.ComputeEdgeFeatures(model_part, edge_index, backend="cupy")
        self.assertEqual(reference.shape, device.shape)
        self.assertTrue(numpy.allclose(reference, device, rtol=0.0, atol=1e-12))


@KratosUnittest.skipUnless(have_cupy, _MISSING)
class TestCuPyBackendStaysOptIn(_ForcedDeviceMixin, KratosUnittest.TestCase):
    """Even with a working GPU and the thresholds at zero, nothing selects
    CuPy unless asked - the guarantee that installing cupy cannot silently
    change anybody's results."""

    def test_AutoStillResolvesToNumpyWithAGpuPresent(self):
        array_backend_utils.SetDefaultArrayBackend(None)
        xp, is_cupy = array_backend_utils.ResolveArrayModule("auto", size_hint=10 ** 9)
        self.assertIs(xp, numpy)
        self.assertFalse(is_cupy)

    def test_TheDefaultArgumentOfAConvertedSiteIsNumpy(self):
        rng = numpy.random.default_rng(6)
        grid = rng.random((2, 5, 5, 5))
        points = rng.random((40, 3))
        box = (numpy.zeros(3), numpy.ones(3))
        self.assertIsInstance(
            grid_bridge.InterpolateGridAtPoints(grid, box, points), numpy.ndarray)


if __name__ == "__main__":
    KratosUnittest.main()
