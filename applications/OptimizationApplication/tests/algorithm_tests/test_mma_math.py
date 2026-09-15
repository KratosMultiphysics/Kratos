"""Validates the core MMA/GCMMA math (algorithms/mma_math.py) against the
closed-form / worked examples published in K. Svanberg, "The Method of
Moving Asymptotes - A New Method for Structural Optimization", Int. J.
Numer. Methods Eng., Vol. 24, 359-373 (1987), Section 6.

This module loads mma_math.py directly by file path (not via the normal
"KratosMultiphysics.OptimizationApplication..." dotted import) so it can run
under plain `python3 -m unittest` without a compiled Kratos build.
"""
import importlib.util
import os
import unittest

import numpy as np


def _load_mma_math():
    module_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "python_scripts", "algorithms", "mma_math.py")
    spec = importlib.util.spec_from_file_location("mma_math", os.path.normpath(module_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mma_math = _load_mma_math()


def _run_mma(objective, constraint, x0, xmin, xmax, num_outer_iter):
    n = x0.size
    m = 1
    x = x0.copy()
    xold1 = x.copy()
    xold2 = x.copy()
    low = None
    upp = None

    a0 = 1.0
    a = np.zeros(m)
    c = np.full(m, 1000.0)
    d = np.full(m, 1.0)

    for outer_iter in range(1, num_outer_iter + 1):
        f0, df0 = objective(x)
        f1_val, df1 = constraint(x)
        f1 = np.array([f1_val])
        df1_mat = df1.reshape(1, -1)

        low, upp = mma_math.update_asymptotes(x, xold1, xold2, xmin, xmax, low, upp, outer_iter)
        alfa, beta = mma_math.compute_move_limits(x, low, upp, xmin, xmax)

        p0, q0 = mma_math.compute_pq(df0, x, low, upp)
        p, q = mma_math.compute_pq(df1_mat, x, low, upp)
        r = mma_math.compute_r(f1, x, p, q, low, upp)
        b = -r

        xmma, _, _, _, _ = mma_math.solve_mma_subproblem(x, alfa, beta, low, upp, p0, q0, p, q, b, a0, a, c, d)

        xold2 = xold1.copy()
        xold1 = x.copy()
        x = xmma.copy()

    return x


def _run_gcmma(objective, constraint, x0, xmin, xmax, num_outer_iter, inner_max_iter=15, tol=1e-6):
    n = x0.size
    m = 1
    x = x0.copy()
    xold1 = x.copy()
    xold2 = x.copy()
    low = None
    upp = None

    a0 = 1.0
    a = np.zeros(m)
    c = np.full(m, 1000.0)
    d = np.full(m, 1.0)

    for outer_iter in range(1, num_outer_iter + 1):
        f0, df0 = objective(x)
        f1_val, df1 = constraint(x)
        f1 = np.array([f1_val])
        df1_mat = df1.reshape(1, -1)

        low, upp = mma_math.update_asymptotes(x, xold1, xold2, xmin, xmax, low, upp, outer_iter)
        alfa, beta = mma_math.compute_move_limits(x, low, upp, xmin, xmax)

        raa0 = mma_math.initial_raa(df0, xmin, xmax)
        raa = mma_math.initial_raa(df1_mat, xmin, xmax)

        xmma = x
        for _ in range(inner_max_iter):
            p0, q0 = mma_math.compute_pq_gcmma(df0, x, xmin, xmax, low, upp, raa0)
            r0 = mma_math.compute_r(f0, x, p0, q0, low, upp)
            p, q = mma_math.compute_pq_gcmma(df1_mat, x, xmin, xmax, low, upp, raa)
            r = mma_math.compute_r(f1, x, p, q, low, upp)
            b = -r

            xmma, _, _, _, _ = mma_math.solve_mma_subproblem(x, alfa, beta, low, upp, p0, q0, p, q, b, a0, a, c, d)

            f0_approx = mma_math.evaluate_approximation(xmma, p0, q0, r0, low, upp)
            f_approx = mma_math.evaluate_approximation(xmma, p, q, r, low, upp)

            f0_real, _ = objective(xmma)
            f1_real_val, _ = constraint(xmma)
            f_real = np.array([f1_real_val])

            if mma_math.check_conservativeness(f0_real, f_real, f0_approx, f_approx, tol):
                break
            raa0, raa = mma_math.update_raa(raa0, raa, f0_real, f_real, f0_approx, f_approx)

        xold2 = xold1.copy()
        xold1 = x.copy()
        x = xmma.copy()

    return x


class TestMMAMathCantileverBeam(unittest.TestCase):
    """Paper Test Problem 1 (Section 6): 5-element cantilever beam with a
    single displacement constraint. Known closed-form optimum (eq. 22):
    x* = [6.016, 5.309, 4.494, 3.502, 2.153], objective = 1.340."""

    C1 = 0.0624
    C2 = 1.0
    COEFFS = np.array([61.0, 37.0, 19.0, 7.0, 1.0])
    EXPECTED_X = np.array([6.016, 5.309, 4.494, 3.502, 2.153])
    EXPECTED_OBJ = 1.340

    @classmethod
    def _objective(cls, x):
        f0 = cls.C1 * np.sum(x)
        df0 = cls.C1 * np.ones_like(x)
        return f0, df0

    @classmethod
    def _constraint(cls, x):
        f1 = np.sum(cls.COEFFS / x ** 3) - cls.C2
        df1 = -3.0 * cls.COEFFS / x ** 4
        return f1, df1

    def setUp(self):
        self.x0 = np.array([5.0, 5.0, 5.0, 5.0, 5.0])
        self.xmin = np.full(5, 0.5)
        self.xmax = np.full(5, 15.0)

    def test_mma_converges_to_paper_optimum(self):
        x = _run_mma(self._objective, self._constraint, self.x0, self.xmin, self.xmax, num_outer_iter=25)

        obj, _ = self._objective(x)
        infeas, _ = self._constraint(x)

        np.testing.assert_allclose(x, self.EXPECTED_X, atol=0.05)
        self.assertAlmostEqual(obj, self.EXPECTED_OBJ, delta=1.0e-3)
        self.assertLessEqual(infeas, 1.0e-3)

    def test_gcmma_converges_to_paper_optimum(self):
        x = _run_gcmma(self._objective, self._constraint, self.x0, self.xmin, self.xmax, num_outer_iter=15)

        obj, _ = self._objective(x)
        infeas, _ = self._constraint(x)

        np.testing.assert_allclose(x, self.EXPECTED_X, atol=0.05)
        self.assertAlmostEqual(obj, self.EXPECTED_OBJ, delta=1.0e-3)
        self.assertLessEqual(infeas, 1.0e-3)


class TestMMAMathTwoBarTruss(unittest.TestCase):
    """Paper Test Problem 3 (Section 6): 2-bar truss with one active stress
    constraint (bar 2's constraint is never active, per the paper, and is
    omitted here). The paper does not publish a single converged optimum for
    this problem -- only an oscillating iterate sequence (Table III) in
    roughly x1 in [1.14, 1.41], x2 in [0.25, 0.50], w in [1.28, 1.53] -- so
    the assertions below are deliberately loose, bounding that envelope
    rather than checking a single closed-form point."""

    C1 = 1.0
    C2 = 0.124

    @classmethod
    def _objective(cls, x):
        x1, x2 = x
        w = cls.C1 * x1 * np.sqrt(1.0 + x2 ** 2)
        dw_dx1 = cls.C1 * np.sqrt(1.0 + x2 ** 2)
        dw_dx2 = cls.C1 * x1 * x2 / np.sqrt(1.0 + x2 ** 2)
        return w, np.array([dw_dx1, dw_dx2])

    @classmethod
    def _constraint(cls, x):
        x1, x2 = x
        s = np.sqrt(1.0 + x2 ** 2)
        term = 8.0 / x1 + 1.0 / (x1 * x2)
        sigma1 = cls.C2 * s * term
        f1 = sigma1 - 1.0

        dterm_dx1 = -8.0 / x1 ** 2 - 1.0 / (x1 ** 2 * x2)
        dterm_dx2 = -1.0 / (x1 * x2 ** 2)
        ds_dx2 = x2 / s

        df1_dx1 = cls.C2 * s * dterm_dx1
        df1_dx2 = cls.C2 * (ds_dx2 * term + s * dterm_dx2)
        return f1, np.array([df1_dx1, df1_dx2])

    def setUp(self):
        self.x0 = np.array([1.5, 0.5])
        self.xmin = np.array([0.2, 0.1])
        self.xmax = np.array([4.0, 1.6])

    def test_mma_stays_within_paper_envelope(self):
        x = _run_mma(self._objective, self._constraint, self.x0, self.xmin, self.xmax, num_outer_iter=25)

        obj, _ = self._objective(x)
        infeas, _ = self._constraint(x)

        self.assertGreaterEqual(x[0], 1.1)
        self.assertLessEqual(x[0], 1.5)
        self.assertGreaterEqual(x[1], 0.2)
        self.assertLessEqual(x[1], 0.55)
        self.assertGreaterEqual(obj, 1.3)
        self.assertLessEqual(obj, 1.6)
        self.assertLessEqual(infeas, 1.0e-3)


if __name__ == "__main__":
    unittest.main()
