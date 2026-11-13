"""Numerical kernels of the Method of Moving Asymptotes (MMA) and its globally convergent variant (GCMMA).

Pure numpy, no Kratos dependency. Implemented from:
    [svanberg-1987-mma]           K. Svanberg, "The method of moving asymptotes - a new method for structural
                                  optimization", Int. J. Numer. Methods Eng. 24, 359-373 (1987).
    [svanberg-2007-mma-gcmma-note] K. Svanberg, "MMA and GCMMA - two methods for nonlinear optimization", KTH note (2007).

Problem solved by the subproblem (standard nonlinear programming, the "a" / "z" extension of the note is not used):
    minimize    f0(x) + sum_i (c_i y_i + 0.5 d_i y_i^2)
    subject to  f_i(x) - y_i <= 0,  xmin <= x <= xmax,  y >= 0.
All arrays are indexed [function, variable]; function 0 is the objective, 1..m are the constraints.
"""
import numpy as np

_TINY = 1e-300


def UpdateAsymptotes(x: np.ndarray,
                     x_prev: 'np.ndarray | None',
                     x_prev2: 'np.ndarray | None',
                     lower_prev: 'np.ndarray | None',
                     upper_prev: 'np.ndarray | None',
                     xmin: np.ndarray,
                     xmax: np.ndarray,
                     asyinit: float = 0.5,
                     asyincr: float = 1.2,
                     asydecr: float = 0.7,
                     asymin: float = 0.01,
                     asymax: float = 10.0) -> 'tuple[np.ndarray, np.ndarray]':
    """Returns the lower and upper asymptotes (l, u) of the current iteration.

    Without two previous points (`x_prev2 is None`) the asymptotes are x -/+ asyinit*(xmax-xmin)
    [svanberg-2007-mma-gcmma-note p.5, eq. 3.11]. Otherwise they are contracted by asydecr where the variable
    oscillates and expanded by asyincr where it moves monotonically (eq. 3.12-3.13), then kept within
    [asymin, asymax]*(xmax-xmin) of x (eq. 3.14).
    """
    span = xmax - xmin
    if x_prev is None or x_prev2 is None or lower_prev is None or upper_prev is None:
        return x - asyinit * span, x + asyinit * span

    sign = (x - x_prev) * (x_prev - x_prev2)
    gamma = np.ones_like(x)
    gamma[sign < 0.0] = asydecr
    gamma[sign > 0.0] = asyincr

    lower = x - gamma * (x_prev - lower_prev)
    upper = x + gamma * (upper_prev - x_prev)

    lower = np.clip(lower, x - asymax * span, x - asymin * span)
    upper = np.clip(upper, x + asymin * span, x + asymax * span)
    return lower, upper


def ComputeMoveLimits(x: np.ndarray,
                      lower: np.ndarray,
                      upper: np.ndarray,
                      xmin: np.ndarray,
                      xmax: np.ndarray,
                      asymptote_margin: float = 0.1,
                      move: float = 0.5) -> 'tuple[np.ndarray, np.ndarray]':
    """Returns the subproblem bounds (alpha, beta) [svanberg-2007-mma-gcmma-note p.4, eq. 3.6-3.7].

    alpha = max(xmin, l + asymptote_margin (x - l), x - move (xmax - xmin)), beta symmetrically with u and xmax.
    asymptote_margin is the public "albefa" setting of AlgorithmMMA.
    """
    span = xmax - xmin
    alpha = np.maximum.reduce([xmin, lower + asymptote_margin * (x - lower), x - move * span])
    beta = np.minimum.reduce([xmax, upper - asymptote_margin * (upper - x), x + move * span])
    return alpha, beta


def ComputeApproximation(x: np.ndarray,
                         f: np.ndarray,
                         df: np.ndarray,
                         lower: np.ndarray,
                         upper: np.ndarray,
                         xmin: np.ndarray,
                         xmax: np.ndarray,
                         regularisation: 'float | np.ndarray' = 1e-5) -> 'tuple[np.ndarray, np.ndarray, np.ndarray]':
    """Returns the MMA approximation coefficients (p, q, r) of f_i(x) ~ sum_j p_ij/(u_j-x_j) + q_ij/(x_j-l_j) + r_i.

    p_ij = (u_j-x_j)^2 (max(g,0) + 0.001|g| + reg_i/(xmax_j-xmin_j)), q_ij = (x_j-l_j)^2 (max(-g,0) + 0.001|g| + ...),
    g = df_ij, and r_i makes the approximation exact at x [svanberg-2007-mma-gcmma-note p.4, eq. 3.3-3.5].
    The 0.001|g| and reg terms keep p, q > 0, hence a strictly convex subproblem. `regularisation` is raa0 (1e-5) for
    MMA and the per-function rho_i for GCMMA (p.7, eq. 4.3-4.4); a scalar or an array of length len(f).

    Args:
        x: Current design, shape (n,).
        f: Function values at x, shape (m+1,).
        df: Gradients at x, shape (m+1, n).
    """
    span = xmax - xmin
    reg = np.broadcast_to(np.asarray(regularisation, dtype=float), f.shape)[:, None] / span[None, :]
    magnitude = 0.001 * np.abs(df)
    p = (upper - x) ** 2 * (np.maximum(df, 0.0) + magnitude + reg)
    q = (x - lower) ** 2 * (np.maximum(-df, 0.0) + magnitude + reg)
    r = f - EvaluateApproximation(x, p, q, np.zeros_like(f), lower, upper)
    return p, q, r


def EvaluateApproximation(x: np.ndarray,
                          p: np.ndarray,
                          q: np.ndarray,
                          r: np.ndarray,
                          lower: np.ndarray,
                          upper: np.ndarray) -> np.ndarray:
    """Evaluates all approximated functions at x, shape (m+1,) [svanberg-2007-mma-gcmma-note p.4, eq. 3.2]."""
    return p @ (1.0 / (upper - x)) + q @ (1.0 / (x - lower)) + r


def _PrimalFromDual(lam: np.ndarray, p: np.ndarray, q: np.ndarray, lower: np.ndarray, upper: np.ndarray,
                    alpha: np.ndarray, beta: np.ndarray):
    """Closed-form minimiser x_j(lam) of the Lagrangian [svanberg-1987-mma p.365, eq. 17-19].

    Returns the clipped x, the unclipped x (to detect the free variables), P_j(lam) and Q_j(lam).
    """
    big_p = p[0] + lam @ p[1:]
    big_q = q[0] + lam @ q[1:]
    root_p = np.sqrt(np.maximum(big_p, _TINY))
    root_q = np.sqrt(np.maximum(big_q, _TINY))
    x_free = (root_p * lower + root_q * upper) / (root_p + root_q)
    return np.clip(x_free, alpha, beta), x_free, big_p, big_q


def SolveSubproblem(p: np.ndarray,
                    q: np.ndarray,
                    b: np.ndarray,
                    lower: np.ndarray,
                    upper: np.ndarray,
                    alpha: np.ndarray,
                    beta: np.ndarray,
                    c: np.ndarray,
                    d: np.ndarray,
                    tolerance: float = 1e-10,
                    max_iter: int = 200) -> 'tuple[np.ndarray, np.ndarray, np.ndarray, float]':
    """Solves the MMA subproblem with the dual method of [svanberg-1987-mma p.364-365, Sec. 4-5].

    minimize g_0(x) + sum_i (c_i y_i + 0.5 d_i y_i^2) s.t. g_i(x) - y_i <= b_i, alpha <= x <= beta, y >= 0, with
    g_i(x) = sum_j p_ij/(u_j-x_j) + q_ij/(x_j-l_j). The Lagrangian is minimised in closed form,
    x_j(lam) clipped to [alpha_j, beta_j] and y_i(lam) = max(0, (lam_i-c_i)/d_i); the concave dual W(lam) is then
    maximised over lam >= 0 by a projected Newton method with Armijo backtracking (the first trial step stops at the
    kink lam_i = c_i of the penalty term). Variables at a bound are exactly alpha/beta.

    Args:
        p, q: Coefficients of all functions, shape (m+1, n); row 0 is the objective.
        b: Right-hand sides of the m constraints (b_i = -r_i).
        c, d: Artificial-variable costs, shape (m,); d must be > 0.
        tolerance: Stop when every projected dual gradient component i is below tolerance*max(|b_i|, |g_i(x)|), i.e.
            relative to the size of each constraint. c must stay commensurate with the constraint scale (Svanberg's
            scaling advice), otherwise the dual stalls at round-off and the residual stays large.

    Returns:
        (x, y, lam, dual_residual) where dual_residual is the largest of these per-constraint relative components
        (1e-6 or less means solved).
    """
    m = b.size
    if np.any(d <= 0.0):
        raise RuntimeError("SolveSubproblem requires d > 0 for every constraint.")
    if m == 0:
        x, _, _, _ = _PrimalFromDual(np.zeros(0), p, q, lower, upper, alpha, beta)
        return x, np.zeros(0), np.zeros(0), 0.0

    p_con, q_con = p[1:], q[1:]

    def Evaluate(lam):
        x, x_free, big_p, big_q = _PrimalFromDual(lam, p, q, lower, upper, alpha, beta)
        excess = np.maximum(lam - c, 0.0)
        g = p_con @ (1.0 / (upper - x)) + q_con @ (1.0 / (x - lower))
        value = np.sum(big_p / (upper - x) + big_q / (x - lower)) - lam @ b - np.sum(excess ** 2 / (2.0 * d))
        gradient = g - b - excess / d
        return value, gradient, x, x_free, big_p, big_q, g

    def ProjectedGradient(lam, gradient):
        return np.where(lam > 0.0, gradient, np.maximum(gradient, 0.0))

    lam = np.zeros(m)
    scale_b = float(np.max(np.abs(b)))

    def Scale(g):
        return max(scale_b, float(np.max(np.abs(g))), _TINY)

    def RelativeResidual(lam, gradient, g):
        # Per-constraint scale, so a small constraint next to a large one is still solved to its own precision.
        constraint_scale = np.maximum(np.maximum(np.abs(b), np.abs(g)), 1e-12 * Scale(g))
        return float(np.max(np.abs(ProjectedGradient(lam, gradient)) / constraint_scale))

    value, gradient, x, x_free, big_p, big_q, g = Evaluate(lam)
    for _ in range(max_iter):
        if RelativeResidual(lam, gradient, g) <= tolerance:
            break

        # Exact Hessian of W on the free variables: -S diag(1/h) S^T - diag(1/d) where y_i > 0 (lam > c, or at the kink
        # lam = c when ascending).
        free = (x_free > alpha) & (x_free < beta)
        slope = p_con / (upper - x) ** 2 - q_con / (x - lower) ** 2
        curvature = 2.0 * big_p / (upper - x) ** 3 + 2.0 * big_q / (x - lower) ** 3
        neg_hessian = (slope[:, free] / curvature[free]) @ slope[:, free].T + np.diag(((lam > c) | ((lam == c) & (gradient > 0.0))) / d)

        at_zero = (lam <= 0.0) & (gradient <= 0.0)
        movable = ~at_zero
        step = np.zeros(m)
        # Levenberg damping relative to the Hessian scale, increased until the Newton step is an ascent direction.
        reg_scale = max(np.max(np.abs(np.diag(neg_hessian))), Scale(g) / (1.0 + np.max(lam)))
        sub_hessian, sub_gradient = neg_hessian[np.ix_(movable, movable)], gradient[movable]
        step[movable] = sub_gradient
        for reg in (1e-14, 1e-10, 1e-6, 1e-2):
            try:
                trial_step = np.linalg.solve(sub_hessian + reg * reg_scale * np.eye(sub_gradient.size), sub_gradient)
            except np.linalg.LinAlgError:
                continue
            if np.all(np.isfinite(trial_step)) and sub_gradient @ trial_step > 0.0:
                step[movable] = trial_step
                break

        # The dual is only piecewise quadratic (kink at lam = c): stop the first trial at the first kink crossed.
        crossing = (lam - c) * (lam + step - c) < 0.0
        kink = int(np.argmin(np.where(crossing, (c - lam) / np.where(step == 0.0, 1.0, step), np.inf))) if crossing.any() else -1
        t = min(1.0, float((c[kink] - lam[kink]) / step[kink])) if kink >= 0 else 1.0
        accepted = False
        for trial in range(200):
            lam_new = np.maximum(lam + t * step, 0.0)
            if trial == 0 and kink >= 0 and t < 1.0:
                lam_new[kink] = c[kink]
            new_value, new_gradient, new_x, new_x_free, new_p, new_q, new_g = Evaluate(lam_new)
            if new_value >= value + 1e-4 * (gradient @ (lam_new - lam)):
                accepted = True
                break
            t *= 0.5
        if not accepted:
            break
        lam, value, gradient, x, x_free, big_p, big_q, g = (
            lam_new, new_value, new_gradient, new_x, new_x_free, new_p, new_q, new_g)

    y = np.maximum(lam - c, 0.0) / d
    return x, y, lam, RelativeResidual(lam, gradient, g)


def KKTResidual(df0: np.ndarray,
                df: np.ndarray,
                f: np.ndarray,
                lam: np.ndarray,
                x: np.ndarray,
                xmin: np.ndarray,
                xmax: np.ndarray,
                bound_tolerance: float = 1e-6) -> float:
    """Returns the KKT residual of the original problem at x: max of stationarity, feasibility, complementarity.

    Stationarity is the inf-norm of the projected gradient of the Lagrangian df0 + lam^T df, where the component of a
    variable within bound_tolerance*(xmax-xmin) of a bound is zeroed if it pushes outwards. Feasibility is max(f, 0) and
    complementarity max |lam_i f_i|.
    """
    gradient = df0 + (lam @ df if lam.size else 0.0)
    tol = bound_tolerance * (xmax - xmin)
    projected = np.where(((x <= xmin + tol) & (gradient > 0.0)) | ((x >= xmax - tol) & (gradient < 0.0)), 0.0, gradient)
    residual = float(np.max(np.abs(projected))) if projected.size else 0.0
    if f.size:
        residual = max(residual, float(np.max(np.maximum(f, 0.0))), float(np.max(np.abs(lam * f))))
    return residual


def InitialRho(df: np.ndarray, xmin: np.ndarray, xmax: np.ndarray) -> np.ndarray:
    """Returns the initial GCMMA parameters rho_i = max(1e-6, 0.1/n sum_j |df_ij| (xmax_j-xmin_j)), shape (m+1,).

    [svanberg-2007-mma-gcmma-note p.7, eq. 4.6]
    """
    return np.maximum(0.1 / df.shape[1] * (np.abs(df) @ (xmax - xmin)), 1e-6)


def RhoScale(x_new: np.ndarray, x: np.ndarray, lower: np.ndarray, upper: np.ndarray, xmin: np.ndarray,
             xmax: np.ndarray) -> float:
    """Returns d(x_new), the coefficient multiplying rho in the GCMMA approximation [svanberg-2007-mma-gcmma-note p.8, eq. 4.7]."""
    return float(np.sum((upper - lower) * (x_new - x) ** 2 / ((upper - x_new) * (x_new - lower) * (xmax - xmin))))


def IsConservative(f_true: np.ndarray, f_approx: np.ndarray, tolerance: float = 1e-6) -> bool:
    """True if f_approx_i >= f_true_i for all i, up to tolerance*max(1, |f_true_i|)
    [svanberg-2007-mma-gcmma-note p.6]."""
    return bool(np.all(f_approx >= f_true - tolerance * np.maximum(1.0, np.abs(f_true))))


def UpdateRho(rho: np.ndarray, f_true: np.ndarray, f_approx: np.ndarray, rho_scale: float,
              tolerance: float = 1e-6) -> np.ndarray:
    """Returns the next rho: min(1.1 (rho + delta), 10 rho) for the non-conservative functions, rho otherwise.

    delta_i = (f_true_i - f_approx_i)/d(x_new) [svanberg-2007-mma-gcmma-note p.8, eq. 4.8-4.9]. A function counts as
    non-conservative with the same tolerance as IsConservative.
    """
    gap = f_true - f_approx
    violated = gap > tolerance * np.maximum(1.0, np.abs(f_true))
    delta = gap / rho_scale if rho_scale > 0.0 else np.full_like(gap, np.inf)
    return np.where(violated, np.minimum(1.1 * (rho + delta), 10.0 * rho), rho)
