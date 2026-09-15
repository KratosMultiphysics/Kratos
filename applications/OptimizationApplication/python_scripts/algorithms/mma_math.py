"""Core numerical routines for the Method of Moving Asymptotes (MMA) and its
globally convergent variant (GCMMA).

This module implements, from first principles, the separable convex
approximation scheme described by K. Svanberg, "The Method of Moving
Asymptotes - A New Method for Structural Optimization", Int. J. Numer.
Methods Eng., Vol. 24, 359-373 (1987), together with the artificial-variable
formulation of Section 5 of that paper, and the general conservativeness
mechanism of GCMMA (Svanberg, 2002).

The convex subproblem generated at every outer iteration is solved with a
primal-dual interior-point Newton method applied directly to the KKT system
of the artificial-variable formulation. This is a standard technique for
convex, separable optimization problems and is not specific to any single
implementation of MMA.

All arrays are plain 1-D numpy arrays: shape (n,) for design-variable
indexed quantities, (m,) for constraint-indexed quantities and (m, n) for
constraint gradient/coefficient matrices. This module has no Kratos
dependency and can be used and tested standalone.
"""

import time

import numpy as np


# ---------------------------------------------------------------------------
# Asymptotes and move limits
# ---------------------------------------------------------------------------

def update_asymptotes(x, xold1, xold2, xmin, xmax, low_prev, upp_prev, outer_iter,
                       asyinit=0.5, asyincr=1.2, asydecr=0.7, asymin=0.01, asymax=10.0):
    """Updates the lower/upper moving asymptotes.

    For the first two outer iterations the fixed-distance rule (paper eq. 9,
    with s0 = asyinit) is used, since xold1/xold2 are not yet available. For
    later iterations the oscillation/monotone heuristic of eq. 11-13 is
    applied: the asymptotes are tightened (asydecr) where the design
    oscillates, and relaxed (asyincr) where it moves monotonically, then
    clipped so they never collapse onto x or run away unboundedly.

    Args:
        x: current design point, shape (n,)
        xold1: design point one outer iteration ago (unused if outer_iter<=2)
        xold2: design point two outer iterations ago (unused if outer_iter<=2)
        xmin, xmax: variable bounds, shape (n,)
        low_prev, upp_prev: asymptotes from the previous outer iteration
            (unused if outer_iter<=2)
        outer_iter: 1-based outer iteration counter
    """
    xmami = xmax - xmin

    if outer_iter <= 2:
        low = x - asyinit * xmami
        upp = x + asyinit * xmami
    else:
        sign_product = (x - xold1) * (xold1 - xold2)
        factor = np.ones_like(x)
        factor[sign_product > 0] = asyincr
        factor[sign_product < 0] = asydecr

        low = x - factor * (xold1 - low_prev)
        upp = x + factor * (upp_prev - xold1)

        low = np.maximum(low, x - asymax * xmami)
        low = np.minimum(low, x - asymin * xmami)
        upp = np.minimum(upp, x + asymax * xmami)
        upp = np.maximum(upp, x + asymin * xmami)

    return low, upp


def compute_move_limits(x, low, upp, xmin, xmax, albefa=0.1, move=0.5):
    """Computes the move limits alfa/beta (paper eq. 8, generalized with a
    move cap and a hard clip to [xmin, xmax])."""
    xmami = xmax - xmin
    alfa = np.maximum(np.maximum(low + albefa * (x - low), x - move * xmami), xmin)
    beta = np.minimum(np.minimum(upp - albefa * (upp - x), x + move * xmami), xmax)
    return alfa, beta


# ---------------------------------------------------------------------------
# Separable approximation coefficients
# ---------------------------------------------------------------------------

def compute_pq(dfdx, x, low, upp):
    """Plain-MMA p, q coefficients (paper eq. 2-4).

    `dfdx` may have shape (n,) (objective gradient) or (m, n) (stacked
    constraint gradients); `x`, `low`, `upp` always have shape (n,). Returns
    p, q with the same shape as `dfdx`.
    """
    pos = np.maximum(dfdx, 0.0)
    neg = np.maximum(-dfdx, 0.0)
    p = pos * (upp - x) ** 2
    q = neg * (x - low) ** 2
    return p, q


def compute_pq_gcmma(dfdx, x, xmin, xmax, low, upp, raa):
    """Conservative (curvature-augmented) p, q coefficients for GCMMA.

    Adds a curvature term proportional to raa/(xmax-xmin) on top of the
    plain-MMA hard positive/negative split, which is the qualitative
    mechanism GCMMA uses to make its convex approximation deliberately more
    conservative than plain MMA's.

    NOTE: the precise mixing constants used by Svanberg's own GCMMA
    formulation (2002) have not been verified against that primary source at
    the time this was written; this uses the straightforward formula above,
    which is mathematically well defined (curvature grows with raa, always
    non-negative) but should be revisited against the primary source before
    being considered final.

    `raa` is either a scalar (for the objective, paired with a 1-D `dfdx`) or
    a vector of shape (m,) (for the constraints, paired with a 2-D `dfdx`).
    """
    xmami = np.maximum(xmax - xmin, 1.0e-5)
    ux1 = upp - x
    xl1 = x - low
    pos = np.maximum(dfdx, 0.0)
    neg = np.maximum(-dfdx, 0.0)

    if np.ndim(dfdx) == 2:
        curvature = np.asarray(raa).reshape(-1, 1) / xmami
    else:
        curvature = raa / xmami

    p = pos * ux1 ** 2 + curvature * ux1 ** 2
    q = neg * xl1 ** 2 + curvature * xl1 ** 2
    return p, q


def compute_r(value, x, p, q, low, upp):
    """Residual constant r (paper eq. 5) such that the approximation matches
    `value` exactly at `x`. `value` is a scalar (with p, q shape (n,)) or a
    vector of shape (m,) (with p, q shape (m, n))."""
    approx_at_x = np.sum(p / (upp - x) + q / (x - low), axis=-1)
    return value - approx_at_x


def evaluate_approximation(x_eval, p, q, r, low, upp):
    """Evaluates the separable approximation r + sum_j(p_j/(upp_j-x_j) +
    q_j/(x_j-low_j)) at an arbitrary point `x_eval`."""
    return r + np.sum(p / (upp - x_eval) + q / (x_eval - low), axis=-1)


# ---------------------------------------------------------------------------
# Primal-dual interior-point subsolve
# ---------------------------------------------------------------------------

_DEFAULT_SUBSOLVE_SETTINGS = {
    "epsilon_init": 1.0,
    "epsilon_min": 1.0e-7,
    "epsilon_reduction_factor": 0.1,
    "residual_tol_factor": 0.9,
    "fraction_to_boundary": 0.99,
    "max_newton_iter": 50,
}


def _compute_residuals(x, y, z, lam, xsi, eta, mu, zeta, s, p0, q0, p, q, b, low, upp,
                        alfa, beta, a0, a, c, d, epsilon):
    ux1 = upp - x
    xl1 = x - low
    plam = p0 + lam @ p
    qlam = q0 + lam @ q
    gvec = p @ (1.0 / ux1) + q @ (1.0 / xl1)

    rex = plam / ux1 ** 2 - qlam / xl1 ** 2 - xsi + eta
    rey = c + d * y - lam - mu
    rez = a0 - a @ lam - zeta
    relam = gvec - a * z - y + s - b
    rexsi = xsi * (x - alfa) - epsilon
    reeta = eta * (beta - x) - epsilon
    remu = mu * y - epsilon
    rezet = zeta * z - epsilon
    res = lam * s - epsilon

    return rex, rey, rez, relam, rexsi, reeta, remu, rezet, res


def _residual_norms(*args):
    parts = _compute_residuals(*args)
    flat = np.concatenate([np.atleast_1d(part) for part in parts])
    return float(np.linalg.norm(flat)), float(np.max(np.abs(flat)))


def _newton_direction(x, y, z, lam, xsi, eta, mu, zeta, s, p0, q0, p, q, b, low, upp,
                       alfa, beta, a0, a, c, d, epsilon):
    ux1 = upp - x
    xl1 = x - low
    plam = p0 + lam @ p
    qlam = q0 + lam @ q
    GG = p / ux1 ** 2 - q / xl1 ** 2

    (rex, rey, rez, relam, rexsi, reeta, remu, rezet, res) = _compute_residuals(
        x, y, z, lam, xsi, eta, mu, zeta, s, p0, q0, p, q, b, low, upp, alfa, beta, a0, a, c, d, epsilon)

    diagx = 2.0 * plam / ux1 ** 3 + 2.0 * qlam / xl1 ** 3 + xsi / (x - alfa) + eta / (beta - x)
    diagy = d + mu / y
    diagz = zeta / z
    diaglam = s / lam

    delx = rex + rexsi / (x - alfa) - reeta / (beta - x)
    dely = rey + remu / y
    delz = rez + rezet / z
    dellam = relam - res / lam

    # The (x, y, z) blocks of the full KKT Jacobian are diagonal (each x_j/y_i only
    # couples to itself and to lam, never to another x_k/y_l), so dx, dy, dz can be
    # eliminated analytically, leaving only an (m x m) system to solve for dlam. This
    # is essential for performance: n is the number of design variables (can be in the
    # thousands for real topology optimization), while m (number of constraints) is
    # typically small -- solving an (m x m) system instead of a dense (n+2m+1) x (n+2m+1)
    # one turns an O(n^3) factorization into O(n*m^2).
    GG_scaled = GG / diagx  # (m, n): each column j divided by diagx_j
    schur = GG_scaled @ GG.T + np.diag(1.0 / diagy) + np.outer(a, a) / diagz + np.diag(diaglam)
    rhs_lam = dellam - GG_scaled @ delx + dely / diagy + a * delz / diagz

    dlam = np.linalg.solve(schur, rhs_lam)
    dz = (-delz + a @ dlam) / diagz
    dy = (-dely + dlam) / diagy
    dx = (-delx - GG.T @ dlam) / diagx

    dxsi = (-rexsi - xsi * dx) / (x - alfa)
    deta = (-reeta + eta * dx) / (beta - x)
    dmu = (-remu - mu * dy) / y
    dzeta = (-rezet - zeta * dz) / z
    ds = (-res - s * dlam) / lam

    return dx, dy, dz, dlam, dxsi, deta, dmu, dzeta, ds


def _max_fraction_to_boundary_step(x, y, z, lam, xsi, eta, mu, zeta, s,
                                    dx, dy, dz, dlam, dxsi, deta, dmu, dzeta, ds,
                                    alfa, beta, tau):
    candidates = [1.0]

    for val, dval in ((y, dy), (lam, dlam), (xsi, dxsi), (eta, deta), (mu, dmu), (s, ds)):
        mask = dval < 0.0
        if np.any(mask):
            candidates.append(np.min(-tau * val[mask] / dval[mask]))

    if dz < 0.0:
        candidates.append(-tau * z / dz)
    if dzeta < 0.0:
        candidates.append(-tau * zeta / dzeta)

    dec = dx < 0.0
    if np.any(dec):
        candidates.append(np.min(-tau * (x[dec] - alfa[dec]) / dx[dec]))
    inc = dx > 0.0
    if np.any(inc):
        candidates.append(np.min(tau * (beta[inc] - x[inc]) / dx[inc]))

    return min(candidates)


def _solve_unconstrained(alfa, beta, p0, q0, low, upp):
    """Closed-form solution of the m=0 (no-constraint) subproblem: each
    variable minimizes its own strictly convex term independently. This is
    equivalent to the case analysis of paper eq. 17-19 (the rational function
    p0_j/(upp_j-x_j) + q0_j/(x_j-low_j) is convex with a single interior
    minimizer, so clipping that minimizer to [alfa_j, beta_j] reproduces the
    three-case check without branching)."""
    sqrt_p0 = np.sqrt(p0)
    sqrt_q0 = np.sqrt(q0)
    denom = sqrt_p0 + sqrt_q0
    with np.errstate(divide="ignore", invalid="ignore"):
        interior = (sqrt_p0 * low + sqrt_q0 * upp) / denom
    x_star = np.where(denom > 0.0, interior, 0.5 * (alfa + beta))
    return np.clip(x_star, alfa, beta)


def solve_mma_subproblem(x0, alfa, beta, low, upp, p0, q0, p, q, b, a0, a, c, d, settings=None, verbose=False):
    """Solves the MMA/GCMMA convex subproblem:

        minimize_{x,y,z}   sum_j(p0_j/(upp_j-x_j) + q0_j/(x_j-low_j))
                            + a0*z + sum_i(c_i*y_i + 0.5*d_i*y_i^2)
        subject to          sum_j(p_ij/(upp_j-x_j) + q_ij/(x_j-low_j)) - a_i*z - y_i <= b_i
                             alfa_j <= x_j <= beta_j,  y_i >= 0,  z >= 0

    via a primal-dual interior-point Newton method on the KKT system of this
    artificial-variable formulation (Svanberg 1987, Section 5).

    If `verbose` is True, prints one line per epsilon level and one line per
    Newton iteration, each with a wall-clock timestamp and elapsed time, so
    progress can be observed on large problems.

    Returns:
        xmma, ymma, zmma, lam, subsolve_residual

    Note: `subsolve_residual` measures how precisely *this iteration's convex
    approximation* was solved (its interior-point continuation always drives
    it down to roughly `subsolve_settings.epsilon_min`, by construction, on
    every outer iteration) -- it is a solver-health diagnostic, not a measure
    of how close the design is to a true optimum of the original problem. For
    that, use `kkt_residual` with the real (unapproximated) gradients/values
    at the accepted design point, together with the `lam` returned here.
    """
    settings = dict(_DEFAULT_SUBSOLVE_SETTINGS, **(settings or {}))
    epsilon = settings["epsilon_init"]
    epsilon_min = settings["epsilon_min"]
    epsilon_reduction_factor = settings["epsilon_reduction_factor"]
    residual_tol_factor = settings["residual_tol_factor"]
    tau = settings["fraction_to_boundary"]
    max_newton_iter = settings["max_newton_iter"]

    m = a.size
    n = x0.size

    if m == 0:
        xmma = _solve_unconstrained(alfa, beta, p0, q0, low, upp)
        return xmma, np.zeros(0), 0.0, np.zeros(0), 0.0

    start_time = time.time()

    def _log(message):
        if verbose:
            elapsed = time.time() - start_time
            print(f"[mma_math] {time.strftime('%H:%M:%S')} (+{elapsed:7.2f}s) {message}", flush=True)

    _log(f"solve_mma_subproblem start: n={n}, m={m}")

    x = 0.5 * (alfa + beta)
    y = np.ones(m)
    z = 1.0
    lam = np.ones(m)
    xsi = np.maximum(1.0 / (x - alfa), 1.0)
    eta = np.maximum(1.0 / (beta - x), 1.0)
    mu = np.maximum(0.5 * c, 1.0)
    zeta = 1.0
    s = np.ones(m)

    epsilon_level = 0
    while epsilon > epsilon_min:
        epsilon_level += 1
        residual_norm, residual_max = _residual_norms(
            x, y, z, lam, xsi, eta, mu, zeta, s, p0, q0, p, q, b, low, upp, alfa, beta, a0, a, c, d, epsilon)
        _log(f"epsilon level {epsilon_level}: epsilon={epsilon:.2e}, residual_max={residual_max:.3e}")

        newton_iter = 0
        while residual_max > residual_tol_factor * epsilon and newton_iter < max_newton_iter:
            newton_iter += 1
            newton_start = time.time()

            dx, dy, dz, dlam, dxsi, deta, dmu, dzeta, ds = _newton_direction(
                x, y, z, lam, xsi, eta, mu, zeta, s, p0, q0, p, q, b, low, upp, alfa, beta, a0, a, c, d, epsilon)

            steg = _max_fraction_to_boundary_step(
                x, y, z, lam, xsi, eta, mu, zeta, s, dx, dy, dz, dlam, dxsi, deta, dmu, dzeta, ds, alfa, beta, tau)

            for _ in range(50):
                x_new = x + steg * dx
                y_new = y + steg * dy
                z_new = z + steg * dz
                lam_new = lam + steg * dlam
                xsi_new = xsi + steg * dxsi
                eta_new = eta + steg * deta
                mu_new = mu + steg * dmu
                zeta_new = zeta + steg * dzeta
                s_new = s + steg * ds

                new_norm, new_max = _residual_norms(
                    x_new, y_new, z_new, lam_new, xsi_new, eta_new, mu_new, zeta_new, s_new,
                    p0, q0, p, q, b, low, upp, alfa, beta, a0, a, c, d, epsilon)
                if new_norm < residual_norm:
                    break
                steg *= 0.5

            x, y, z, lam = x_new, y_new, z_new, lam_new
            xsi, eta, mu, zeta, s = xsi_new, eta_new, mu_new, zeta_new, s_new
            residual_norm, residual_max = new_norm, new_max

            _log(f"  newton iter {newton_iter}: step={steg:.3e}, residual_max={residual_max:.3e}, "
                 f"iter_time={time.time() - newton_start:.3f}s")

        epsilon *= epsilon_reduction_factor

    _, kkt_max = _residual_norms(
        x, y, z, lam, xsi, eta, mu, zeta, s, p0, q0, p, q, b, low, upp, alfa, beta, a0, a, c, d, 0.0)

    _log(f"solve_mma_subproblem done: subsolve_residual={kkt_max:.3e}, total_time={time.time() - start_time:.3f}s")

    return x, y, z, lam, kkt_max


def kkt_residual(x, df0dx, dfdx, f_value, xmin, xmax, lam):
    """True KKT residual of the *original*, unapproximated problem at design
    point `x`, using the real objective/constraint gradients and values and
    the Lagrange multipliers `lam` returned by the last-solved MMA/GCMMA
    subproblem.

    Combines three standard KKT ingredients for `minimize f0(x) subject to
    f_i(x) <= 0, xmin <= x <= xmax`:
      - stationarity: the box-constrained projected-gradient residual of the
        Lagrangian `grad f0(x) + lam @ grad f(x)` (zero iff x is a KKT point):
        for each variable strictly interior to [xmin, xmax] the raw gradient
        component must vanish; at the lower bound only a negative component
        counts against it (a positive one just means the bound is correctly
        binding), and symmetrically at the upper bound. This has no arbitrary
        step-size parameter -- unlike a naive "take a unit step and clip to
        the box" residual, which saturates at (and gets stuck reporting) a
        value of order (xmax-xmin) whenever the raw gradient is larger than
        the box width, regardless of how close x actually is to optimal;
      - feasibility: how much any f_i(x) exceeds 0;
      - complementary slackness: how far lam_i * f_i(x) is from 0.

    Deliberately NOT rescaled by the gradient's own norm (e.g. dividing
    df0dx by ||df0dx||_inf): since the max-abs component of any nonzero
    vector divided by its own max-abs is always exactly 1, that would make
    this residual pinned at ~1 for any nonzero gradient and only ever drop
    to 0 at *exact* convergence -- a binary signal, not a genuinely shrinking
    one. So this residual is intentionally reported on the same raw scale as
    the gradients passed in (i.e. the responses' own scaling), and is not
    directly comparable in magnitude to a differently-normalized residual
    from another algorithm.

    Unlike the subsolve's own interior-point residual (see
    `solve_mma_subproblem`), this is not artificially bounded by a fixed
    solver tolerance -- it genuinely shrinks as the outer MMA/GCMMA
    iterations converge to a stationary point of the real problem.
    """
    grad_lagrangian = df0dx + lam @ dfdx if lam.size else df0dx
    at_lower = x <= xmin
    at_upper = x >= xmax
    projected_grad = np.where(at_lower, np.minimum(grad_lagrangian, 0.0), grad_lagrangian)
    projected_grad = np.where(at_upper, np.maximum(projected_grad, 0.0), projected_grad)
    stationarity = np.max(np.abs(projected_grad)) if x.size else 0.0

    if f_value.size:
        feasibility = np.max(np.maximum(f_value, 0.0))
        complementarity = np.max(np.abs(lam * f_value))
    else:
        feasibility = 0.0
        complementarity = 0.0

    return max(stationarity, feasibility, complementarity)


def initial_raa(dfdx, xmin, xmax, factor=0.1, floor=1.0e-6):
    """A scale-aware starting curvature parameter for GCMMA.

    The curvature contributed by `raa` to the p, q coefficients is
    raa/(xmax-xmin) (see `compute_pq_gcmma`); for that contribution to be a
    meaningful fraction of the approximation from the first inner iteration
    (rather than needing dozens of doublings to become numerically
    significant), `raa` itself should scale with gradient magnitude times
    variable range -- this is a basic dimensional argument, not a tuned
    constant. `dfdx` is 1-D (objective, scalar result) or 2-D (constraints,
    one result per row).
    """
    xmami = np.maximum(xmax - xmin, 1.0e-5)
    if np.ndim(dfdx) == 2:
        scale = factor * np.mean(np.abs(dfdx) * xmami, axis=-1)
    else:
        scale = factor * np.mean(np.abs(dfdx) * xmami)
    return np.maximum(scale, floor)


# ---------------------------------------------------------------------------
# GCMMA conservativeness bookkeeping
# ---------------------------------------------------------------------------

def check_conservativeness(f0_real, f_real, f0_approx, f_approx, tol=1.0e-6):
    """True iff the approximation over-estimates (or matches, within `tol`)
    both the real objective and every real constraint value at the trial
    point -- the GCMMA acceptance condition."""
    if f0_approx < f0_real - tol:
        return False
    if np.size(f_real) and np.any(f_approx < f_real - tol):
        return False
    return True


def update_raa(raa0, raa, f0_real, f_real, f0_approx, f_approx, growth_factor=2.0, raa_max=1.0e5):
    """Grows the curvature parameters for the objective (raa0) and for every
    constraint whose approximation under-estimated the real value, via
    simple geometric backoff. This is a deliberate simplification of
    Svanberg's own (2002) minimal-sufficient-increase rule, guaranteed to
    terminate since unbounded curvature growth eventually dominates any
    bounded nonlinearity on a compact box."""
    raa0_new = min(raa0 * growth_factor, raa_max) if f0_approx < f0_real else raa0

    raa_new = np.array(raa, dtype=float, copy=True)
    if np.size(raa_new):
        violated = f_approx < f_real
        raa_new[violated] = np.minimum(raa_new[violated] * growth_factor, raa_max)

    return raa0_new, raa_new
