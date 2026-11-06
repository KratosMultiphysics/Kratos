"""The array-backend selector: numpy stays the default, CuPy is opt-in.

The first class runs everywhere, including the torch-free CI - selection is
pure Python and must behave with no cupy installed. The second only runs
where a CUDA device actually answers.
"""

import os

import numpy

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils

try:
    import cupy
    have_cupy = array_backend_utils.IsCuPyAvailable()
except ImportError:
    have_cupy = False


class TestArrayBackendSelection(KratosUnittest.TestCase):
    """Selection logic, which must work with or without cupy present."""

    def setUp(self):
        array_backend_utils.SetDefaultArrayBackend(None)
        self._environment = os.environ.pop(
            array_backend_utils._ENVIRONMENT_VARIABLE, None)

    def tearDown(self):
        array_backend_utils.SetDefaultArrayBackend(None)
        os.environ.pop(array_backend_utils._ENVIRONMENT_VARIABLE, None)
        if self._environment is not None:
            os.environ[array_backend_utils._ENVIRONMENT_VARIABLE] = self._environment

    def test_AutoIsNumpy(self):
        """The whole point of the default: nobody gets CuPy by accident."""
        xp, is_cupy = array_backend_utils.ResolveArrayModule("auto")
        self.assertIs(xp, numpy)
        self.assertFalse(is_cupy)

    def test_NumpyIsAlwaysHonoured(self):
        xp, is_cupy = array_backend_utils.ResolveArrayModule("numpy", size_hint=10 ** 9)
        self.assertIs(xp, numpy)
        self.assertFalse(is_cupy)

    def test_UnknownBackendRaises(self):
        with self.assertRaises(ValueError):
            array_backend_utils.ResolveArrayModule("gpu")
        with self.assertRaises(ValueError):
            array_backend_utils.SetDefaultArrayBackend("gpu")

    def test_BelowThresholdFallsBackToNumpy(self):
        """A small problem never goes to the device, even when asked."""
        xp, is_cupy = array_backend_utils.ResolveArrayModule("cupy", size_hint=10)
        self.assertIs(xp, numpy)
        self.assertFalse(is_cupy)

    def test_PerSiteThresholdIsHonoured(self):
        xp, _ = array_backend_utils.ResolveArrayModule(
            "cupy", size_hint=200000, threshold=5000000)
        self.assertIs(xp, numpy)

    def test_EnvironmentVariableSetsTheDefault(self):
        os.environ[array_backend_utils._ENVIRONMENT_VARIABLE] = "numpy"
        self.assertEqual(array_backend_utils.GetDefaultArrayBackend(), "numpy")
        self.assertIs(array_backend_utils.ResolveArrayModule("auto")[0], numpy)

    def test_AnUnparseableEnvironmentValueIsIgnored(self):
        os.environ[array_backend_utils._ENVIRONMENT_VARIABLE] = "nonsense"
        self.assertEqual(array_backend_utils.GetDefaultArrayBackend(), "auto")

    def test_SetDefaultOverridesTheEnvironment(self):
        os.environ[array_backend_utils._ENVIRONMENT_VARIABLE] = "cupy"
        array_backend_utils.SetDefaultArrayBackend("numpy")
        self.assertEqual(array_backend_utils.GetDefaultArrayBackend(), "numpy")

    def test_AvailabilityIsABoolNeverAnException(self):
        """Called on every converted path, so it must not raise on a CPU box."""
        self.assertIsInstance(array_backend_utils.IsCuPyAvailable(), bool)

    def test_ToHostAcceptsNumpy(self):
        values = numpy.arange(6.0).reshape(2, 3)
        self.assertTrue(numpy.array_equal(array_backend_utils.ToHost(values), values))


@KratosUnittest.skipUnless(have_cupy, "Missing required python module: cupy (with a CUDA device).")
class TestArrayBackendCuPy(KratosUnittest.TestCase):
    """The device paths, where a GPU is actually present."""

    def test_CupyIsSelectedAboveTheThreshold(self):
        xp, is_cupy = array_backend_utils.ResolveArrayModule("cupy", size_hint=10 ** 8)
        self.assertTrue(is_cupy)
        self.assertIsNot(xp, numpy)

    def test_ToDeviceAndBackRoundTrips(self):
        values = numpy.arange(12.0).reshape(3, 4)
        xp, _ = array_backend_utils.ResolveArrayModule("cupy", size_hint=10 ** 8)
        on_device = array_backend_utils.ToDevice(values, xp)
        self.assertFalse(isinstance(on_device, numpy.ndarray))
        back = array_backend_utils.ToHost(on_device)
        self.assertIsInstance(back, numpy.ndarray)
        self.assertTrue(numpy.array_equal(back, values))

    def test_ToDeviceWithNumpyIsANoOp(self):
        values = numpy.arange(4.0)
        self.assertIsInstance(array_backend_utils.ToDevice(values, numpy), numpy.ndarray)

    def test_AsTorchTensorKeepsACupyArrayOnTheDevice(self):
        """The reason the helper exists: no device-to-host-to-device bounce."""
        try:
            import torch
        except ImportError:
            self.skipTest("Missing required python module: torch.")
        xp, _ = array_backend_utils.ResolveArrayModule("cupy", size_hint=10 ** 8)
        values = xp.arange(12.0).reshape(3, 4)
        tensor = array_backend_utils.AsTorchTensor(values)
        self.assertEqual(tensor.device.type, "cuda")
        self.assertTrue(numpy.allclose(tensor.cpu().numpy(),
                                       array_backend_utils.ToHost(values)))

    def test_AsTorchTensorFromNumpyStaysOnTheHost(self):
        try:
            import torch
        except ImportError:
            self.skipTest("Missing required python module: torch.")
        tensor = array_backend_utils.AsTorchTensor(numpy.arange(6.0))
        self.assertEqual(tensor.device.type, "cpu")


if __name__ == "__main__":
    KratosUnittest.main()
