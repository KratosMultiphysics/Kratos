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
        with self.assertRaisesRegex(Exception, 'Unknown array backend'):
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


class TestCupyShimsWithMocks(KratosUnittest.TestCase):
    """Verifies the cupyx shim signatures without a GPU by stubbing the cupyx modules.

    Guards the review blocker: cupyx gmres is keyword-only "rtol" (no "tol"), and
    cupyx solve_triangular accepts "check_finite".
    """
    def _InstallFakeCupyx(self, sparse_linalg, linalg=None):
        import sys, types
        cupyx = types.ModuleType("cupyx")
        scipy = types.ModuleType("cupyx.scipy")
        sparse = types.ModuleType("cupyx.scipy.sparse")
        sparse_linalg_mod = types.ModuleType("cupyx.scipy.sparse.linalg")
        sparse_linalg_mod.gmres = sparse_linalg
        linalg_mod = types.ModuleType("cupyx.scipy.linalg")
        if linalg is not None:
            linalg_mod.solve_triangular = linalg
        saved = {k: sys.modules.get(k) for k in ["cupyx", "cupyx.scipy", "cupyx.scipy.sparse", "cupyx.scipy.sparse.linalg", "cupyx.scipy.linalg"]}
        sys.modules["cupyx"] = cupyx
        sys.modules["cupyx.scipy"] = scipy
        sys.modules["cupyx.scipy.sparse"] = sparse
        sys.modules["cupyx.scipy.sparse.linalg"] = sparse_linalg_mod
        sys.modules["cupyx.scipy.linalg"] = linalg_mod
        return saved

    def _RestoreModules(self, saved):
        import sys
        for key, mod in saved.items():
            if mod is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = mod

    def test_gmres_cupy_uses_rtol(self):
        from KratosMultiphysics.CoSimulationApplication.utilities.array_backend import ArrayBackend
        calls = {}
        def fake_gmres(op, rhs, **kwargs):
            calls.update(kwargs)
            return ("ok", 0)
        saved = self._InstallFakeCupyx(fake_gmres)
        try:
            backend = ArrayBackend(xp=None, name="cupy")
            result, info = backend.Gmres(object(), np.ones(2), atol=1e-14, rtol=1e-5)
            self.assertEqual(result, "ok")
            self.assertIn("rtol", calls)
            self.assertNotIn("tol", calls)
            self.assertAlmostEqual(calls["rtol"], 1e-5)
        finally:
            self._RestoreModules(saved)

    def test_gmres_cupy_falls_back_to_tol_for_old_cupyx(self):
        from KratosMultiphysics.CoSimulationApplication.utilities.array_backend import ArrayBackend
        calls = []
        def fake_gmres(op, rhs, **kwargs):
            calls.append(dict(kwargs))
            if "rtol" in kwargs:
                raise TypeError("unexpected keyword 'rtol'")
            return ("ok-old", 0)
        saved = self._InstallFakeCupyx(fake_gmres)
        try:
            backend = ArrayBackend(xp=None, name="cupy")
            result, _ = backend.Gmres(object(), np.ones(2), atol=1e-14, rtol=1e-5)
            self.assertEqual(result, "ok-old")
            self.assertEqual(len(calls), 2)
            self.assertIn("rtol", calls[0])
            self.assertIn("tol", calls[1])
        finally:
            self._RestoreModules(saved)

    def test_solve_triangular_cupy_passes_check_finite(self):
        from KratosMultiphysics.CoSimulationApplication.utilities.array_backend import ArrayBackend
        calls = {}
        def fake_solve_triangular(a, b, **kwargs):
            calls.update(kwargs)
            return "ok-tri"
        saved = self._InstallFakeCupyx(lambda *a, **k: (None, None), linalg=fake_solve_triangular)
        try:
            backend = ArrayBackend(xp=None, name="cupy")
            result = backend.SolveTriangular(np.eye(2), np.ones(2))
            self.assertEqual(result, "ok-tri")
            self.assertIn("check_finite", calls)
            self.assertFalse(calls["check_finite"])
        finally:
            self._RestoreModules(saved)


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

    def test_iqnils_cupy_matches_numpy(self):
        self._SkipIfNoCupy()
        expected = self._RunSimpleAccelerator("iqnils", "numpy")
        actual = self._RunSimpleAccelerator("iqnils", "cupy")
        for exp, act in zip(expected, actual):
            np.testing.assert_allclose(act, exp, atol=1e-10)

    def test_anderson_cupy_matches_numpy(self):
        self._SkipIfNoCupy()
        expected = self._RunSimpleAccelerator("anderson", "numpy")
        actual = self._RunSimpleAccelerator("anderson", "cupy")
        for exp, act in zip(expected, actual):
            np.testing.assert_allclose(act, exp, atol=1e-10)

    def test_block_mvqn_cupy_matches_numpy(self):
        self._SkipIfNoCupy()
        expected = self._RunBlockAccelerator("block_mvqn", "numpy")
        actual = self._RunBlockAccelerator("block_mvqn", "cupy")
        for exp, act in zip(expected, actual):
            np.testing.assert_allclose(act, exp, atol=1e-10)

    def test_aitken_identical_residuals_does_not_raise(self):
        # regression: two bit-identical residuals give a zero denominator; the old
        # np.float64 division produced nan/inf and continued, it must not raise
        # ZeroDivisionError (see review on float() placement in aitken.py)
        settings = KM.Parameters('{ "type" : "aitken", "backend" : "numpy" }')
        acc = CreateConvergenceAccelerator(settings)
        acc.InitializeSolutionStep()
        r = np.array([1.0, 2.0, 3.0])
        x = np.array([0.1, 0.2, 0.3])
        acc.UpdateSolution(r, x)
        with np.errstate(invalid="ignore"): # 0/0 -> nan is the expected numpy semantic here
            delta = acc.UpdateSolution(r.copy(), x.copy()) # r_diff == 0
        self.assertIsInstance(delta, np.ndarray)
        self.assertEqual(delta.shape, r.shape)
        self.assertTrue(np.all(np.isnan(delta))) # nan propagates, run continues (no exception)
        acc.FinalizeSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
