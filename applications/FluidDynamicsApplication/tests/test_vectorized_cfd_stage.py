import numpy as np
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication.cfd_utils as cfd_utils_module
import KratosMultiphysics.FluidDynamicsApplication.python_pool as python_pool_module

# Configure the backend for unit-testing the stage without instantiating a simulation.
# If another test module already configured it, reuse that configuration.
try:
    cfd_utils_module.configure("OpenMP", "float64")
    _reset_backend_on_exit = True
except RuntimeError:
    _reset_backend_on_exit = False

import KratosMultiphysics.FluidDynamicsApplication.vectorized_cfd_stage as vcs_module

# The stage module routes array calls through a module-level backend alias normally set inside
# VectorizedCFDStage.__init__; mirror it here for the __init__-free tests below.
vcs_module.xp = cfd_utils_module.xp


def tearDownModule():
    if _reset_backend_on_exit:
        cfd_utils_module.reset()


class TestVectorizedCFDStageHelpers(UnitTest.TestCase):
    """Unit tests for VectorizedCFDStage helper methods, exercised on lightweight stand-ins."""

    def _make_stage(self):
        # Instantiate without running __init__ (which builds a full simulation backend):
        # the helpers under test only read the attributes each test sets explicitly.
        return object.__new__(vcs_module.VectorizedCFDStage)

    # ------------------------------------------------------------------
    # _CheckRejectionLimit
    # ------------------------------------------------------------------
    def test__CheckRejectionLimit_allows_up_to_limit(self):
        stage = self._make_stage()
        stage.max_consecutive_rejections = 3
        stage.cfl = 0.5
        for n in range(1, 4):
            stage._CheckRejectionLimit(n, "step 1")  # must not raise

    def test__CheckRejectionLimit_raises_beyond_limit(self):
        stage = self._make_stage()
        stage.max_consecutive_rejections = 3
        stage.cfl = 0.5
        with self.assertRaises(RuntimeError):
            stage._CheckRejectionLimit(4, "step 2")

    # ------------------------------------------------------------------
    # _ComputeDeltaTime
    # ------------------------------------------------------------------
    def _reference_dt(self, stage, rho_conv, max_dt_fourier=7.0, user_dt=10.0):
        eps = np.finfo(rho_conv.dtype).eps
        dt_cfl = stage.cfl / (rho_conv[np.argmax(rho_conv)] + eps)
        if stage.is_startup:
            dt_cfl *= 0.5
        return float(min(min(dt_cfl, max_dt_fourier), user_dt))

    def test__ComputeDeltaTime_no_startup_reference(self):
        stage = self._make_stage()
        stage.cfl = 0.5
        stage.dt = 10.0
        stage.dt_fourier = 7.0
        stage.is_startup = False
        rho_conv = np.array([1.0, 4.0, 0.5])
        max_cfl_el_id, dt = stage._ComputeDeltaTime(rho_conv)
        self.assertEqual(int(max_cfl_el_id), 2)  # 1-based index of the max-CFL element
        self.assertAlmostEqual(dt, self._reference_dt(stage, rho_conv), delta=1e-12)

    def test__ComputeDeltaTime_startup_halves_dt_cfl(self):
        # Regression: 'st_cfl *= 0.5' used to raise NameError whenever startup_time > 0.
        stage = self._make_stage()
        stage.cfl = 0.5
        stage.dt = 10.0
        stage.dt_fourier = 7.0
        rho_conv = np.array([1.0, 4.0, 0.5])

        stage.is_startup = False
        _, dt_no_startup = stage._ComputeDeltaTime(rho_conv)
        stage.is_startup = True
        _, dt_startup = stage._ComputeDeltaTime(rho_conv)

        self.assertAlmostEqual(dt_startup, 0.5 * dt_no_startup, delta=1e-12)
        self.assertAlmostEqual(dt_startup, self._reference_dt(stage, rho_conv), delta=1e-12)

    def test__ComputeDeltaTime_clamped_by_fourier_bound(self):
        stage = self._make_stage()
        stage.cfl = 0.5
        stage.dt = 10.0
        stage.dt_fourier = 7.0
        stage.is_startup = False
        _, dt = stage._ComputeDeltaTime(np.zeros(4))
        self.assertAlmostEqual(dt, 7.0, delta=1e-12)

    def test__ComputeDeltaTime_clamped_by_user_dt(self):
        stage = self._make_stage()
        stage.cfl = 0.5
        stage.dt = 2.0
        stage.dt_fourier = 9.0
        stage.is_startup = False
        _, dt = stage._ComputeDeltaTime(np.zeros(4))
        self.assertAlmostEqual(dt, 2.0, delta=1e-12)

    def test__ComputeDeltaTime_rejects_invalid_input(self):
        stage = self._make_stage()
        stage.cfl = 0.5
        stage.dt = 10.0
        stage.dt_fourier = 7.0
        stage.is_startup = False
        with self.assertRaises(ValueError):
            stage._ComputeDeltaTime(np.zeros(0))
        with self.assertRaises(ValueError):
            stage._ComputeDeltaTime(np.array([1.0, -0.5]))

    # ------------------------------------------------------------------
    # ComputePressureProjection
    # ------------------------------------------------------------------
    def _make_projection_stage(self):
        rng = np.random.RandomState(1234)
        stage = self._make_stage()
        stage.dim = 2
        stage.n_in_el = 3
        stage.nnodes = 4
        stage.DN = rng.randn(3, 3, 2)
        stage.N = np.ones(3) / 3.0
        stage.elemental_volumes = np.array([0.4, 0.5, 0.6])
        stage.connectivity = np.array([[0, 1, 2], [0, 2, 3], [1, 2, 3]], dtype=np.int32)
        stage.Minv = 1.0 / (1.0 + np.arange(4))   # nodal scalar form (nnodes,)
        stage.cfd_utils = cfd_utils_module.CFDUtils()
        stage.v_el = np.empty(stage.DN.shape)
        stage.v_el_mean_scratch = np.empty((stage.DN.shape[0], stage.dim))
        stage.p_el = np.empty((stage.DN.shape[0], stage.n_in_el))
        stage.pool = python_pool_module.BufferPool(
            sizes=[3 * 3 * 3, 3 * 3 * 2, 3 * 3 * 2], xp=np, dtype=np.float64)
        return stage

    def test_ComputePressureProjection_rejects_wrong_out_shape(self):
        # Regression: the error path used to raise NameError ('pi_press' undefined) instead of ValueError
        stage = self._make_projection_stage()
        with self.assertRaises(ValueError):
            stage.ComputePressureProjection(np.zeros((3, 3)), out=np.zeros((4, 3)))

    def test_ComputePressureProjection_reference_and_pool_release(self):
        stage = self._make_projection_stage()
        rng = np.random.RandomState(42)
        p_elemental = rng.randn(3, 3)

        out = np.zeros((stage.nnodes, stage.dim))
        result = stage.ComputePressureProjection(p_elemental, out=out)

        # Reference: out[e,I,k] = N[I] * sum_J DN[e,J,k] pel[e,J], scaled by volume, assembled, scaled by Minv
        proj_el = np.einsum("I,eJk,eJ->eIk", stage.N, stage.DN, p_elemental)
        proj_el *= stage.elemental_volumes[:, None, None]
        nodal = np.zeros((stage.nnodes, stage.dim))
        np.add.at(nodal, stage.connectivity, proj_el)
        ref = nodal * stage.Minv[:, None]  # Minv nodal, broadcast over dim

        np.testing.assert_allclose(np.asarray(result), ref, rtol=1e-12, atol=1e-14)

        # Every pool buffer acquired inside must be released again (no dangling in_use flags)
        self.assertFalse(any(stage.pool.Status()))

    def test_ComputePressureProjection_requires_out(self):
        stage = self._make_projection_stage()
        with self.assertRaises(ValueError):
            stage.ComputePressureProjection(np.zeros((3, 3)))

    # ------------------------------------------------------------------
    # ComputeVelocityProjection / ComputeDivergenceProjection
    # ------------------------------------------------------------------
    def test_ComputeVelocityProjection_reference_and_pool_release(self):
        stage = self._make_projection_stage()
        rng = np.random.RandomState(7)
        v_elemental = rng.randn(3, 3, 2)

        out = np.zeros((stage.nnodes, stage.dim))
        result = stage.ComputeVelocityProjection(v_elemental, out)

        grad_v = np.einsum("eik,eil->ekl", v_elemental, stage.DN)
        a_grad = np.einsum("eil,ekl->eik", v_elemental, grad_v)
        M_e = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]) / 12.0
        conv = np.einsum("IJ,eJl->eIl", M_e, a_grad)
        conv = conv * stage.elemental_volumes[:, None, None]
        nodal = np.zeros((stage.nnodes, stage.dim))
        np.add.at(nodal, stage.connectivity, conv)
        ref = nodal * stage.Minv[:, None]  # Minv nodal, broadcast over dim

        np.testing.assert_allclose(np.asarray(result), ref, rtol=1e-10, atol=1e-12)
        self.assertFalse(any(stage.pool.Status()))

    def test_ComputeDivergenceProjection_reference_and_pool_release(self):
        stage = self._make_projection_stage()
        rng = np.random.RandomState(11)
        v_elemental = rng.randn(3, 3, 2)

        out = np.zeros(stage.nnodes)
        result = stage.ComputeDivergenceProjection(v_elemental, out=out)

        div = np.einsum("eIk,eIk->e", stage.DN, v_elemental)
        nodal = np.zeros(stage.nnodes)
        np.add.at(nodal, stage.connectivity, div[:, None] * stage.N[None, :] * stage.elemental_volumes[:, None])
        ref = nodal * stage.Minv  # Minv nodal

        np.testing.assert_allclose(np.asarray(result), ref, rtol=1e-10, atol=1e-12)
        self.assertFalse(any(stage.pool.Status()))

    # ------------------------------------------------------------------
    # ComputeElementalConvectiveOperatorSpectralRadius (feeds the M1 substep flow)
    # ------------------------------------------------------------------
    def test_ComputeElementalConvectiveOperatorSpectralRadius_reference_and_pool_release(self):
        stage = self._make_projection_stage()
        rng = np.random.RandomState(21)
        v = rng.randn(stage.nnodes, stage.dim)

        out = stage.ComputeElementalConvectiveOperatorSpectralRadius(v)

        v_el = np.take(v, stage.connectivity, axis=0)
        v_mean = v_el.sum(axis=1) / stage.n_in_el
        proj = np.einsum("ek,eIk->eI", v_mean, stage.DN)
        ref = np.max(np.abs(proj), axis=1)

        np.testing.assert_allclose(np.asarray(out), ref, rtol=1e-10, atol=1e-12)
        self.assertFalse(any(stage.pool.Status()))

    # ------------------------------------------------------------------
    # Normals: one-time normalization (ComputeUnitNormalsFromGlobal) + BC physics
    # ------------------------------------------------------------------
    def test_ComputeUnitNormalsFromGlobal_reference(self):
        stage = self._make_stage()
        ids = np.array([1, 2, 4])
        normals = np.zeros((5, 2))
        normals[1] = np.array([0.0, 2.0])    # non-unit -> normalized
        normals[2] = np.array([-1.0, 0.0])   # already unit
        # row 4 stays zero -> reported as zero-normal and left untouched by clipping
        unit, zero_ids = stage.ComputeUnitNormalsFromGlobal(normals, ids)
        np.testing.assert_allclose(unit, np.array([[0.0, 1.0], [-1.0, 0.0], [0.0, 0.0]]),
                                   rtol=0.0, atol=1e-14)
        np.testing.assert_array_equal(zero_ids, np.array([4]))

    def test_ComputeUnitNormalsFromGlobal_zero_threshold(self):
        # Norm below 1e-15 counts as zero-normal (but is still divided by the 1e-15 clip);
        # norm at/above 1e-15 normalizes to a true unit vector.
        stage = self._make_stage()
        normals = np.array([[1e-16, 0.0], [2e-15, 0.0]])
        unit, zero_ids = stage.ComputeUnitNormalsFromGlobal(normals, np.array([0, 1]))
        np.testing.assert_array_equal(zero_ids, np.array([0]))
        np.testing.assert_allclose(unit[1], np.array([1.0, 0.0]), rtol=0.0, atol=1e-14)

    def test_ApplyOutletBackflowCorrection_removes_only_inflow(self):
        stage = self._make_stage()
        stage.outlet_node_ids = np.array([1, 2, 3, 4])

        normals = np.zeros((5, 2))
        normals[1:4] = np.array([1.0, 0.0])  # unit +x normals; node 4 (outlet) keeps a zero normal
        v = np.array([[5.0, 5.0],   # node 0: not an outlet node -> untouched
                      [2.0, 1.0],   # outflow -> untouched
                      [-3.0, 0.5],  # inflow -> normal component removed
                      [0.0, 4.0],   # tangential -> untouched
                      [7.0, 7.0]])  # zero-normal outlet node -> velocity forced to zero

        # Production-like setup: normals normalized once, zero-normal nodes detected once
        stage.outlet_unit_normals, stage.outlet_zero_normal_indices = \
            stage.ComputeUnitNormalsFromGlobal(normals, stage.outlet_node_ids)

        stage.ApplyOutletBackflowCorrection(v)

        ref = np.array([[5.0, 5.0],
                        [2.0, 1.0],
                        [0.0, 0.5],
                        [0.0, 4.0],
                        [0.0, 0.0]])
        np.testing.assert_allclose(v, ref, rtol=0.0, atol=1e-12)

        # Second call is a no-op (idempotent: outflow/tangential rows untouched again)
        stage.ApplyOutletBackflowCorrection(v)
        np.testing.assert_allclose(v, ref, rtol=0.0, atol=1e-12)

        # The stored boundary normals are unit-length where defined (normalized once, untouched after)
        np.testing.assert_allclose(
            np.linalg.norm(np.asarray(stage.outlet_unit_normals)[:3], axis=1),
            np.ones(3), rtol=1e-12, atol=1e-12)

    def test_ApplyVelocitySlipConditions_removes_normal_component(self):
        stage = self._make_stage()
        stage.slip_node_ids = np.array([1, 2])

        normals = np.zeros((4, 2))
        normals[1] = np.array([0.0, 2.0])   # non-unit normal normalized by the one-time setup
        normals[2] = np.array([1.0, 0.0])
        v = np.array([[9.0, 9.0],
                      [1.0, -3.0],  # strip normal (+y) component -> [1, 0]
                      [2.0, 4.0]])  # strip normal (+x) component -> [0, 4]
        v = np.vstack([v, np.array([[0.0, 0.0]])])  # pad to nnodes=4

        stage.slip_unit_normals, stage.slip_zero_normal_indices = \
            stage.ComputeUnitNormalsFromGlobal(normals, stage.slip_node_ids)

        stage.ApplyVelocitySlipConditions(v)

        ref = np.array([[9.0, 9.0],
                        [1.0, 0.0],
                        [0.0, 4.0],
                        [0.0, 0.0]])
        np.testing.assert_allclose(v, ref, rtol=0.0, atol=1e-12)

    # ------------------------------------------------------------------
    # Buffer pool protocol relied on by the substep loop
    # ------------------------------------------------------------------
    def test_pool_double_get_raises_guarding_the_substep_protocol(self):
        # The substep loop keeps pool buffer 2 acquired (rho_conv) through its last read
        # (ComputeTau) and releases it right before SolveStep1 re-acquires it. Any premature
        # re-acquisition aliases rho_conv's memory and must be caught loudly, never silently.
        pool = python_pool_module.BufferPool(sizes=[8, 8, 8], xp=np, dtype=np.float64)
        rho_conv = pool.Get(2, (4,))
        with self.assertRaises(RuntimeError):
            pool.Get(2, (4,))
        pool.Release(rho_conv)
        # After the release the buffer is available again
        reacquired = pool.Get(2, (4,))
        self.assertTrue(pool.in_use[2])
        pool.Release(reacquired)


if __name__ == '__main__':
    UnitTest.main()
