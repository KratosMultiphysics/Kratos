import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.utilities.array_backend import GetArrayBackend, HostArrayBoundary, IsCupyAvailable
from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator

import numpy as np


class TestArrayBackendResolution(KratosUnittest.TestCase):

    def test_numpy_backend(self):
        backend = GetArrayBackend("numpy")
        self.assertEqual(backend.name, "numpy")
        self.assertIs(backend.xp, np)

    def test_unknown_backend_raises(self):
        with self.assertRaises(Exception):
            GetArrayBackend("not_a_backend")

    def test_cupy_backend_falls_back_to_numpy_if_unavailable(self):
        if IsCupyAvailable():
            self.skipTest("cupy is available on this machine, cannot test the fallback path")
        # explicitly requesting "cupy" without cupy being usable must not raise, it falls back to numpy
        backend = GetArrayBackend("cupy", context="UnitTest")
        self.assertEqual(backend.name, "numpy")

    def test_auto_backend_resolves_to_numpy_without_cupy(self):
        if IsCupyAvailable():
            self.skipTest("cupy is available on this machine, cannot test the numpy-only fallback")
        backend = GetArrayBackend("auto")
        self.assertEqual(backend.name, "numpy")

    def test_cupy_backend_is_usable_if_available(self):
        if not IsCupyAvailable():
            self.skipTest("cupy is not available on this machine")
        backend = GetArrayBackend("cupy")
        self.assertEqual(backend.name, "cupy")


class _DummyAccelerator:
    """Minimal stand-in exposing the "self.backend" attribute the decorator relies on."""
    def __init__(self, backend):
        self.backend = backend

    @HostArrayBoundary
    def UpdateSolution(self, r, x):
        return r * 2.0 + x * 0.0 # depends only on "r", to keep the check simple


class TestHostArrayBoundaryDecorator(KratosUnittest.TestCase):

    def test_passthrough_with_numpy_backend(self):
        dummy = _DummyAccelerator(GetArrayBackend("numpy"))
        r = np.array([1.0, 2.0, 3.0])
        x = np.array([0.1, 0.2, 0.3])
        result = dummy.UpdateSolution(r, x)
        self.assertIsInstance(result, np.ndarray)
        np.testing.assert_array_equal(result, r * 2.0)

    def test_upload_download_roundtrip_with_cupy_backend(self):
        if not IsCupyAvailable():
            self.skipTest("cupy is not available on this machine")
        dummy = _DummyAccelerator(GetArrayBackend("cupy"))
        r = np.array([1.0, 2.0, 3.0])
        x = np.array([0.1, 0.2, 0.3])
        result = dummy.UpdateSolution(r, x)
        # the decorator must download the result back to host numpy
        self.assertIsInstance(result, np.ndarray)
        np.testing.assert_array_equal(result, r * 2.0)


class TestConvergenceAcceleratorsWithBackendSetting(KratosUnittest.TestCase):
    """Directly exercises the numerical convergence accelerators with the "backend" setting.
    These accelerators have no other direct unit tests (elsewhere only the wrapper is tested,
    with the accelerator itself mocked), so this also guards the numpy/cupy conversion.
    """
    def setUp(self):
        self.r = np.array([1.0, 2.0, 3.0])
        self.x = np.array([0.1, 0.2, 0.3])
        self.r2 = np.array([1.1, 1.9, 3.2])
        self.x2 = np.array([0.12, 0.19, 0.31])
        self.y = np.array([0.5, 0.6, 0.7])
        self.y2 = np.array([0.52, 0.55, 0.75])
        self.y_residual = np.array([0.05, 0.02, -0.01])
        self.y_residual2 = np.array([0.04, 0.03, -0.02])

    def _RunSimpleAccelerator(self, acc_type, backend):
        settings = KM.Parameters("""{
            "type"    : "%s",
            "backend" : "%s"
        }""" % (acc_type, backend))
        acc = CreateConvergenceAccelerator(settings)
        self.assertEqual(acc.backend.name, backend if backend != "cupy" or IsCupyAvailable() else "numpy")
        acc.InitializeSolutionStep()
        delta_1 = acc.UpdateSolution(self.r, self.x)
        self.assertIsInstance(delta_1, np.ndarray)
        delta_2 = acc.UpdateSolution(self.r2, self.x2)
        self.assertIsInstance(delta_2, np.ndarray)
        acc.FinalizeSolutionStep()
        return delta_1, delta_2

    def _RunBlockAccelerator(self, acc_type, backend):
        settings = KM.Parameters("""{
            "type"    : "%s",
            "backend" : "%s",
            "solver_sequence" : [
                { "solver" : "a", "data_name" : "data_1", "coupled_data_name" : "data_2" },
                { "solver" : "b", "data_name" : "data_2", "coupled_data_name" : "data_1" }
            ]
        }""" % (acc_type, backend))
        acc = CreateConvergenceAccelerator(settings)
        acc.InitializeSolutionStep()
        delta_1 = acc.UpdateSolution(self.r, self.x, self.y, "data_1", self.y_residual)
        self.assertIsInstance(delta_1, np.ndarray)
        # a distinct "y" is used here so the coupled-solver difference matrix isn't singular,
        # which exercises the actual QR/pinv_product code path (incl. the scipy/cupyx shim)
        delta_2 = acc.UpdateSolution(self.r2, self.x2, self.y2, "data_1", self.y_residual2)
        self.assertIsInstance(delta_2, np.ndarray)
        acc.FinalizeSolutionStep()
        return delta_1, delta_2

    def _SkipIfNoCupy(self):
        if not IsCupyAvailable():
            self.skipTest("cupy is not available on this machine")

    # --- default backend (regression: the new "backend" setting must not change existing behavior) ---

    def test_default_backend_is_numpy(self):
        settings = KM.Parameters('{ "type" : "mvqn" }')
        acc = CreateConvergenceAccelerator(settings)
        self.assertEqual(acc.backend.name, "numpy")

    # --- numpy backend (must always work, no optional dependency involved) ---

    def test_mvqn_numpy(self):          self._RunSimpleAccelerator("mvqn", "numpy")
    def test_aitken_numpy(self):        self._RunSimpleAccelerator("aitken", "numpy")
    def test_iqnils_numpy(self):        self._RunSimpleAccelerator("iqnils", "numpy")
    def test_anderson_numpy(self):      self._RunSimpleAccelerator("anderson", "numpy")
    def test_block_mvqn_numpy(self):    self._RunBlockAccelerator("block_mvqn", "numpy")
    def test_block_ibqnls_numpy(self):  self._RunBlockAccelerator("block_ibqnls", "numpy")

    # --- cupy backend (skipped when cupy / a CUDA device is not available) ---

    def test_mvqn_cupy_matches_numpy(self):
        self._SkipIfNoCupy()
        expected = self._RunSimpleAccelerator("mvqn", "numpy")
        actual = self._RunSimpleAccelerator("mvqn", "cupy")
        for exp, act in zip(expected, actual):
            np.testing.assert_allclose(act, exp, atol=1e-10)

    def test_aitken_cupy_matches_numpy(self):
        self._SkipIfNoCupy()
        expected = self._RunSimpleAccelerator("aitken", "numpy")
        actual = self._RunSimpleAccelerator("aitken", "cupy")
        for exp, act in zip(expected, actual):
            np.testing.assert_allclose(act, exp, atol=1e-10)

    def test_block_ibqnls_cupy_matches_numpy(self):
        self._SkipIfNoCupy()
        expected = self._RunBlockAccelerator("block_ibqnls", "numpy")
        actual = self._RunBlockAccelerator("block_ibqnls", "cupy")
        for exp, act in zip(expected, actual):
            np.testing.assert_allclose(act, exp, atol=1e-10)


if __name__ == '__main__':
    KratosUnittest.main()
