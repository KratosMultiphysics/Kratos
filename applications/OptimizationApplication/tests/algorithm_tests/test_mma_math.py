import numpy as np

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.algorithms import mma_math


def _RunMMA(function, x0, xmin, xmax, constraint_count, iterations, move=0.5, c=1000.0, d=1.0, gcmma=False, inner_max=30):
    """Plain MMA / GCMMA outer loop on a function returning (values, gradients). Returns x, last kkt, last y."""
    x = x0.copy()
    x_prev = x_prev2 = lower = upper = None
    c_vec, d_vec = np.full(constraint_count, c), np.full(constraint_count, d)
    kkt, y = np.inf, np.zeros(constraint_count)
    for _ in range(iterations):
        f, df = function(x)
        lower, upper = mma_math.UpdateAsymptotes(x, x_prev, x_prev2, lower, upper, xmin, xmax)
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, move)
        rho = mma_math.InitialRho(df, xmin, xmax)
        for _ in range(inner_max if gcmma else 1):
            p, q, r = mma_math.ComputeApproximation(x, f, df, lower, upper, xmin, xmax, rho if gcmma else 1e-5)
            x_new, y, lam, _ = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, c_vec, d_vec)
            if not gcmma:
                break
            f_approx = mma_math.EvaluateApproximation(x_new, p, q, r, lower, upper)
            f_true = function(x_new)[0]
            if mma_math.IsConservative(f_true, f_approx):
                break
            rho = mma_math.UpdateRho(rho, f_true, f_approx, mma_math.RhoScale(x_new, x, lower, upper, xmin, xmax))
        kkt = mma_math.KKTResidual(df[0], df[1:], f[1:], lam, x, xmin, xmax)
        x_prev2, x_prev, x = x_prev, x, x_new
    return x, kkt, y


_CANTILEVER_C = np.array([61.0, 37.0, 19.0, 7.0, 1.0])


def _Cantilever(x):
    """Test problem 1 of [svanberg-1987-mma p.366]: min 0.0624 sum(x) s.t. sum(C_j / x_j^3) <= 1."""
    f = np.array([0.0624 * x.sum(), (_CANTILEVER_C / x ** 3).sum() - 1.0])
    df = np.vstack([0.0624 * np.ones(5), -3.0 * _CANTILEVER_C / x ** 4])
    return f, df


def _TwoBarTruss(x):
    """Test problem 3 [svanberg-1987-mma p.372-373]: min x1 sqrt(1+x2^2) s.t. 0.124 sqrt(1+x2^2) (8/x1 +- 1/(x1 x2)) <= 1."""
    s = np.sqrt(1.0 + x[1] ** 2)
    f = [x[0] * s]
    df = [[s, x[0] * x[1] / s]]
    for sigma in (1.0, -1.0):
        h = s * (8.0 * x[1] + sigma) / (x[0] * x[1])
        f.append(0.124 * h - 1.0)
        df.append([-0.124 * h / x[0], 0.124 * ((8.0 * x[1] + sigma) / (s * x[0]) - sigma * s / (x[0] * x[1] ** 2))])
    return np.array(f), np.array(df)


class TestMMAMathApproximation(kratos_unittest.TestCase):
    def test_ComputeApproximationHandValues(self):
        # one variable: x = 1, l = 0, u = 3, xmin = 0, xmax = 2, f = 4, df = 2
        x, lower, upper = np.array([1.0]), np.array([0.0]), np.array([3.0])
        xmin, xmax = np.array([0.0]), np.array([2.0])
        p, q, r = mma_math.ComputeApproximation(x, np.array([4.0]), np.array([[2.0]]), lower, upper, xmin, xmax, 1e-5)
        self.assertAlmostEqual(p[0, 0], 4.0 * (2.0 + 0.002 + 0.5e-5), 12)
        self.assertAlmostEqual(q[0, 0], 1.0 * (0.002 + 0.5e-5), 12)
        self.assertAlmostEqual(r[0], 4.0 - p[0, 0] / 2.0 - q[0, 0] / 1.0, 12)

    def test_ComputeApproximationIsFirstOrderExact(self):
        rng = np.random.default_rng(3)
        n, m = 4, 2
        x = rng.uniform(0.2, 0.8, n)
        xmin, xmax = np.zeros(n), np.ones(n)
        lower, upper = x - rng.uniform(0.5, 1.5, n), x + rng.uniform(0.5, 1.5, n)
        f, df = rng.normal(size=m + 1), rng.normal(size=(m + 1, n))
        for regularisation in (1e-5, np.array([0.1, 1.0, 10.0])):
            p, q, r = mma_math.ComputeApproximation(x, f, df, lower, upper, xmin, xmax, regularisation)
            self.assertTrue(np.all(p > 0.0) and np.all(q > 0.0))
            self.assertVectorAlmostEqual(mma_math.EvaluateApproximation(x, p, q, r, lower, upper), f, 12)
            slope = p / (upper - x) ** 2 - q / (x - lower) ** 2
            self.assertLess(np.max(np.abs(slope - df)), 1e-12)

    def test_ZeroGradientDoesNotDrift(self):
        # With df = 0 only the regularisation is left; it must keep x (not the asymptote mid-point).
        n = 3
        x = np.array([0.2, 0.5, 0.9])
        xmin, xmax = np.zeros(n), np.ones(n)
        lower, upper = x - 0.5, x + 0.5
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax)
        p, q, r = mma_math.ComputeApproximation(x, np.array([1.0, -1.0]), np.zeros((2, n)), lower, upper, xmin, xmax)
        x_new, y, lam, _ = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, np.array([1000.0]), np.array([1.0]))
        self.assertVectorAlmostEqual(x_new, x, 12)
        self.assertVectorAlmostEqual(y, [0.0])
        self.assertVectorAlmostEqual(lam, [0.0])

    def test_ComputeMoveLimits(self):
        x, xmin, xmax = np.array([2.0, 5.0]), np.zeros(2), np.full(2, 10.0)
        lower, upper = np.array([-3.0, 0.0]), np.array([7.0, 6.0])
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.5)
        # alpha_0 = max(0, -2.5, -3) = 0 (variable bound), alpha_1 = max(0, 0.5, 0) = 0.5 (asymptote)
        self.assertVectorAlmostEqual(alpha, [0.0, 0.5])
        # beta_0 = min(10, 6.5, 7) = 6.5 (asymptote), beta_1 = min(10, 5.9, 10) = 5.9
        self.assertVectorAlmostEqual(beta, [6.5, 5.9])
        # a small move limit dominates
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.05)
        self.assertVectorAlmostEqual(alpha, [1.5, 4.5])
        self.assertVectorAlmostEqual(beta, [2.5, 5.5])


class TestMMAMathAsymptotes(kratos_unittest.TestCase):
    def setUp(self):
        self.xmin, self.xmax = np.zeros(3), np.ones(3)

    def test_InitialAsymptotes(self):
        x = np.array([0.2, 0.5, 0.9])
        lower, upper = mma_math.UpdateAsymptotes(x, None, None, None, None, self.xmin, self.xmax, asyinit=0.5)
        self.assertVectorAlmostEqual(lower, x - 0.5)
        self.assertVectorAlmostEqual(upper, x + 0.5)

    def test_SecondIterationUsesInitialRule(self):
        # x_prev known but x_prev2 is None (e.g. after a restart): must not fail and must use the initial rule.
        x = np.array([0.2, 0.5, 0.9])
        lower, upper = mma_math.UpdateAsymptotes(x, x - 0.1, None, x - 1.0, x + 1.0, self.xmin, self.xmax, asyinit=0.3)
        self.assertVectorAlmostEqual(lower, x - 0.3)
        self.assertVectorAlmostEqual(upper, x + 0.3)

    def test_IncreaseAndDecreaseBranches(self):
        x = np.full(3, 0.6)
        x_prev = np.full(3, 0.4)
        x_prev2 = np.array([0.2, 0.6, 0.4])  # monotone, oscillating, zero previous step
        lower_prev, upper_prev = x_prev - 0.5, x_prev + 0.5
        lower, upper = mma_math.UpdateAsymptotes(x, x_prev, x_prev2, lower_prev, upper_prev, self.xmin, self.xmax)
        # gamma = 1.2, 0.7, 1: l = x - gamma 0.5, u = x + gamma 0.5
        self.assertVectorAlmostEqual(lower, [0.0, 0.25, 0.1])
        self.assertVectorAlmostEqual(upper, [1.2, 0.95, 1.1])

    def test_AsymptoteDistanceIsClamped(self):
        x, x_prev, x_prev2 = np.full(3, 0.6), np.full(3, 0.4), np.full(3, 0.2)
        # a huge previous distance is capped at asymax (10) and a tiny one is raised to asymin (0.01)
        lower, upper = mma_math.UpdateAsymptotes(x, x_prev, x_prev2, x_prev - 100.0, x_prev + 0.001, self.xmin, self.xmax)
        self.assertVectorAlmostEqual(lower, x - 10.0)
        self.assertVectorAlmostEqual(upper, x + 0.01)


class TestMMAMathSubproblem(kratos_unittest.TestCase):
    def test_TwoConstraintsAgainstReference(self):
        x = np.array([1.0, 2.0, 3.0])
        xmin, xmax = np.zeros(3), np.full(3, 5.0)
        f = np.array([1.0, 0.5, -0.2])
        df = np.array([[1.0, -2.0, 0.5], [-1.0, 0.5, 1.0], [0.3, -1.0, -1.0]])
        lower, upper = x - 2.5, x + 2.5
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.3)
        p, q, r = mma_math.ComputeApproximation(x, f, df, lower, upper, xmin, xmax)
        x_new, y, lam, residual = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, np.full(2, 1000.0), np.ones(2))

        # reference from scipy trust-constr (xtol 1e-14) on the primal subproblem with y
        self.assertVectorAlmostEqual(x_new, [1.63866527, 2.41598607, 2.73294422], 6)
        self.assertVectorAlmostEqual(y, [0.0, 0.0])
        self.assertLess(residual, 1e-10)
        # both approximated constraints are active, with positive multipliers
        approx = mma_math.EvaluateApproximation(x_new, p, q, r, lower, upper)
        self.assertVectorAlmostEqual(approx[1:], [0.0, 0.0], 9)
        self.assertTrue(np.all(lam > 1.0))
        # stationarity of the Lagrangian for the free variables
        gradient = p[0] / (upper - x_new) ** 2 - q[0] / (x_new - lower) ** 2
        for i in range(2):
            gradient += lam[i] * (p[i + 1] / (upper - x_new) ** 2 - q[i + 1] / (x_new - lower) ** 2)
        self.assertVectorAlmostEqual(gradient, np.zeros(3), 9)

    def test_ScaleInvariance(self):
        # Scaling all functions by s scales p, q, r and leaves x, lam unchanged; the stop test must be relative.
        x = np.array([1.0, 2.0, 3.0])
        xmin, xmax = np.zeros(3), np.full(3, 5.0)
        f = np.array([1.0, 0.5, -0.2])
        df = np.array([[1.0, -2.0, 0.5], [-1.0, 0.5, 1.0], [0.3, -1.0, -1.0]])
        lower, upper = x - 2.5, x + 2.5
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.3)
        for scale in (1e-6, 1e4):
            p, q, r = mma_math.ComputeApproximation(x, scale * f, scale * df, lower, upper, xmin, xmax, 1e-5 * scale)
            x_new, y, lam, residual = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, np.full(2, 1000.0), np.ones(2))
            self.assertVectorAlmostEqual(x_new, [1.63866527, 2.41598607, 2.73294422], 6)
            self.assertLess(residual, 1e-6)
            approx = mma_math.EvaluateApproximation(x_new, p, q, r, lower, upper)
            self.assertLess(np.max(np.abs(approx[1:])), 1e-6 * scale)

    def test_MixedConstraintScales(self):
        # Constraints 1e6 apart in scale must both be solved to their own precision (c large, so y stays 0).
        x = np.array([1.0, 2.0, 3.0])
        xmin, xmax = np.zeros(3), np.full(3, 5.0)
        scales = np.array([1.0, 1e3, 1e-3])
        f = scales * np.array([1.0, 0.5, -0.2])
        df = scales[:, None] * np.array([[1.0, -2.0, 0.5], [-1.0, 0.5, 1.0], [0.3, -1.0, -1.0]])
        lower, upper = x - 2.5, x + 2.5
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.3)
        p, q, r = mma_math.ComputeApproximation(x, f, df, lower, upper, xmin, xmax, 1e-5 * scales)
        x_new, y, lam, residual = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, np.full(2, 1e9), np.ones(2))
        self.assertVectorAlmostEqual(y, [0.0, 0.0])
        self.assertVectorAlmostEqual(x_new, [1.63866527, 2.41598607, 2.73294422], 6)
        self.assertLess(residual, 1e-6)
        approx = mma_math.EvaluateApproximation(x_new, p, q, r, lower, upper)
        self.assertLess(abs(approx[2]), 1e-8 * scales[2])

    def test_VariablesAtBoundAreExact(self):
        # min x1 + x2 with a violated constraint 1 - x1 <= 0: x1 -> beta_1, x2 -> alpha_2 exactly.
        x = np.array([0.5, 0.5])
        xmin, xmax = np.full(2, 0.1), np.full(2, 10.0)
        lower, upper = x - 5.0, x + 5.0
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.01)
        f, df = np.array([1.0, 0.5]), np.array([[1.0, 1.0], [-1.0, 0.0]])
        p, q, r = mma_math.ComputeApproximation(x, f, df, lower, upper, xmin, xmax)
        x_new, y, lam, _ = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, np.array([1000.0]), np.array([1.0]))
        self.assertEqual(x_new[0], beta[0])
        self.assertEqual(x_new[1], alpha[1])

    def test_InfeasibleProblemUsesArtificialVariable(self):
        # constraint 2 - x1 <= 0 with x1 <= 1 cannot be satisfied: y > 0 and x1 sits at its bound.
        x = np.array([0.5, 0.5])
        xmin, xmax = np.zeros(2), np.ones(2)
        lower, upper = x - 1.0, x + 1.0
        alpha, beta = mma_math.ComputeMoveLimits(x, lower, upper, xmin, xmax, 0.1, 0.5)
        f, df = np.array([1.0, 1.5]), np.array([[1.0, 1.0], [-1.0, 0.0]])
        p, q, r = mma_math.ComputeApproximation(x, f, df, lower, upper, xmin, xmax)
        x_new, y, lam, residual = mma_math.SolveSubproblem(p, q, -r[1:], lower, upper, alpha, beta, np.array([1000.0]), np.array([1.0]))
        self.assertEqual(x_new[0], beta[0])
        self.assertGreater(y[0], 0.4)
        approx = mma_math.EvaluateApproximation(x_new, p, q, r, lower, upper)
        self.assertAlmostEqual(approx[1], y[0], 8)
        self.assertLess(residual, 1e-8)

    def test_InvalidArtificialCoefficient(self):
        p, q = np.ones((2, 1)), np.ones((2, 1))
        with self.assertRaisesRegex(RuntimeError, "d > 0"):
            mma_math.SolveSubproblem(p, q, np.zeros(1), -np.ones(1), 2.0 * np.ones(1), np.zeros(1), np.ones(1), np.ones(1), np.zeros(1))


class TestMMAMathKKT(kratos_unittest.TestCase):
    def setUp(self):
        self.xmin, self.xmax = np.full(2, 0.1), np.full(2, 10.0)

    def test_BoundActiveVariableIsZeroed(self):
        # min x1 + x2 s.t. 1 - x1 <= 0 at the optimum [1, 0.1]: x2 at its lower bound with a positive gradient.
        df0, df, f = np.array([1.0, 1.0]), np.array([[-1.0, 0.0]]), np.array([0.0])
        self.assertAlmostEqual(mma_math.KKTResidual(df0, df, f, np.array([1.0]), np.array([1.0, 0.1]), self.xmin, self.xmax), 0.0, 12)

    def test_InteriorGradientCounts(self):
        df0, df, f = np.array([1.0, 1.0]), np.array([[-1.0, 0.0]]), np.array([0.0])
        self.assertAlmostEqual(mma_math.KKTResidual(df0, df, f, np.array([1.0]), np.array([1.0, 0.5]), self.xmin, self.xmax), 1.0, 12)
        # at the lower bound but with a gradient pushing inwards (negative): not zeroed
        self.assertAlmostEqual(mma_math.KKTResidual(-df0, df, f, np.array([0.0]), np.array([1.0, 0.1]), self.xmin, self.xmax), 1.0, 12)

    def test_FeasibilityAndComplementarity(self):
        df0, df = np.zeros(2), np.array([[0.0, 0.0]])
        x = np.array([5.0, 5.0])
        self.assertAlmostEqual(mma_math.KKTResidual(df0, df, np.array([0.3]), np.array([0.0]), x, self.xmin, self.xmax), 0.3, 12)
        self.assertAlmostEqual(mma_math.KKTResidual(df0, df, np.array([-0.5]), np.array([2.0]), x, self.xmin, self.xmax), 1.0, 12)
        self.assertAlmostEqual(mma_math.KKTResidual(df0, np.zeros((0, 2)), np.zeros(0), np.zeros(0), x, self.xmin, self.xmax), 0.0, 12)


class TestMMAMathGCMMA(kratos_unittest.TestCase):
    def test_InitialRho(self):
        df = np.array([[1.0, -3.0], [0.0, 0.0]])
        rho = mma_math.InitialRho(df, np.zeros(2), np.array([2.0, 4.0]))
        self.assertVectorAlmostEqual(rho, [0.1 / 2.0 * (1.0 * 2.0 + 3.0 * 4.0), 1e-6])

    def test_RhoScale(self):
        value = mma_math.RhoScale(np.array([2.0]), np.array([1.0]), np.array([0.0]), np.array([4.0]), np.array([0.0]), np.array([4.0]))
        self.assertAlmostEqual(value, 4.0 * 1.0 / (2.0 * 2.0 * 4.0), 14)

    def test_IsConservative(self):
        f_true = np.array([1.0, 5.0, -1.0])
        self.assertTrue(mma_math.IsConservative(f_true, np.array([1.2, 5.0, -1.0])))
        self.assertTrue(mma_math.IsConservative(f_true, np.array([1.0, 5.0 - 4e-6, -1.0]), 1e-6))   # within 1e-6 max(1, 5)
        self.assertFalse(mma_math.IsConservative(f_true, np.array([1.0, 5.0 - 6e-6, -1.0]), 1e-6))
        # the tolerance is relative to max(1, |f|): tiny functions use the absolute 1e-6
        self.assertFalse(mma_math.IsConservative(np.array([1e-3]), np.array([1e-3 - 2e-6]), 1e-6))

    def test_UpdateRho(self):
        rho = np.array([0.5, 0.5, 0.5])
        f_true, f_approx = np.array([1.0, 5.0, -1.0]), np.array([1.2, 4.0, -1.0])
        # only function 1 is non-conservative: delta = 1/0.25 = 4, min(1.1 (0.5 + 4), 10 0.5) = 4.95
        new_rho = mma_math.UpdateRho(rho, f_true, f_approx, 0.25)
        self.assertVectorAlmostEqual(new_rho, [0.5, 4.95, 0.5], 12)
        # a huge gap is limited to a tenfold increase
        self.assertVectorAlmostEqual(mma_math.UpdateRho(rho, f_true, f_approx, 0.01), [0.5, 5.0, 0.5], 12)
        # d(x) = 0 cannot divide: tenfold increase
        self.assertVectorAlmostEqual(mma_math.UpdateRho(rho, f_true, f_approx, 0.0), [0.5, 5.0, 0.5], 12)


class TestMMAMathCantileverBeam(kratos_unittest.TestCase):
    """[svanberg-1987-mma p.366]: optimum x = (6.016, 5.309, 4.494, 3.502, 2.153), weight 1.340 (reference to 8 digits from SLSQP)."""
    _optimum = np.array([6.01601588, 5.30917385, 4.49432956, 3.50147498, 2.15266535])

    def test_MMA(self):
        x, kkt, y = _RunMMA(_Cantilever, np.full(5, 5.0), np.ones(5), np.full(5, 20.0), 1, 40)
        self.assertVectorAlmostEqual(x, self._optimum, 5)
        self.assertAlmostEqual(_Cantilever(x)[0][0], 1.3399563606, 8)
        self.assertLess(kkt, 1e-6)
        self.assertVectorAlmostEqual(y, [0.0])

    def test_GCMMAWithWideBounds(self):
        # GCMMA converges with a wide box and a large move limit.
        xmin, xmax = np.full(5, 0.1), np.full(5, 50.0)
        x, kkt, _ = _RunMMA(_Cantilever, np.full(5, 5.0), xmin, xmax, 1, 40, move=0.2, gcmma=True)
        self.assertVectorAlmostEqual(x, self._optimum, 5)
        self.assertLess(kkt, 1e-6)


class TestMMAMathTwoBarTruss(kratos_unittest.TestCase):
    _optimum = np.array([1.41163114, 0.37707243])   # SLSQP, ftol 1e-15
    _bounds = (np.array([0.2, 0.1]), np.array([4.0, 1.6]))

    def test_GradientsAgainstFiniteDifferences(self):
        x = np.array([1.3, 0.45])
        _, df = _TwoBarTruss(x)
        for j in range(2):
            step = np.zeros(2)
            step[j] = 1e-6
            fd = (_TwoBarTruss(x + step)[0] - _TwoBarTruss(x - step)[0]) / 2e-6
            self.assertVectorAlmostEqual(df[:, j], fd, 7)

    def test_MMA(self):
        x, kkt, _ = _RunMMA(_TwoBarTruss, np.array([1.5, 0.5]), *self._bounds, 2, 30)
        self.assertVectorAlmostEqual(x, self._optimum, 5)
        self.assertAlmostEqual(_TwoBarTruss(x)[0][0], 1.50865242, 7)
        self.assertLess(kkt, 1e-6)

    def test_GCMMA(self):
        x, kkt, _ = _RunMMA(_TwoBarTruss, np.array([1.5, 0.5]), *self._bounds, 2, 30, gcmma=True)
        self.assertVectorAlmostEqual(x, self._optimum, 5)
        self.assertLess(kkt, 1e-6)


class TestMMAMathBoundActiveProblem(kratos_unittest.TestCase):
    def test_KKTGoesToZero(self):
        # min x1 + x2 s.t. 1 - x1 <= 0, x in [0.1, 10]: optimum [1, 0.1] with x2 at its bound.
        def function(x):
            return np.array([x.sum(), 1.0 - x[0]]), np.array([[1.0, 1.0], [-1.0, 0.0]])
        x, kkt, y = _RunMMA(function, np.array([5.0, 5.0]), np.full(2, 0.1), np.full(2, 10.0), 1, 30)
        self.assertVectorAlmostEqual(x, [1.0, 0.1], 8)
        self.assertLess(kkt, 1e-8)
        self.assertVectorAlmostEqual(y, [0.0])


if __name__ == "__main__":
    kratos_unittest.main()
