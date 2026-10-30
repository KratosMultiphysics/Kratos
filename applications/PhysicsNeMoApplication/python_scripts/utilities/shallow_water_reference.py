"""Synthetic linear shallow-water reference data for grid surrogates.

A pure-numpy finite-difference integrator of the LINEARIZED shallow-water
equations on a doubly-periodic lat-lon-like grid,

    dh/dt = -H0 (du/dx + dv/dy)
    du/dt = -g dh/dx
    dv/dt = -g dh/dy

with centered differences (numpy.roll) and second-order midpoint (RK2)
time stepping under a CFL guard. It exists to feed GraphCast-style
grid-to-grid surrogates with physical step pairs when
ShallowWaterApplication is not compiled - the state layout (3, H, W) =
(height anomaly, u, v) matches GraphCastNet's (C, H, W) grid contract
directly (see the Graph Neural Networks documentation page).

numpy-only and eagerly importable: no torch, no physicsnemo.
"""

import numpy

GRAVITY = 9.81


def _Derivative(state, g, depth, dx, dy):
    """Time derivative of the (3, H, W) state (periodic, centered)."""
    height, u, v = state

    def ddx(field):
        return (numpy.roll(field, -1, axis=1) - numpy.roll(field, 1, axis=1)) / (2.0 * dx)

    def ddy(field):
        return (numpy.roll(field, -1, axis=0) - numpy.roll(field, 1, axis=0)) / (2.0 * dy)

    return numpy.stack([
        -depth * (ddx(u) + ddy(v)),
        -g * ddx(height),
        -g * ddy(height),
    ])


def Step(state, dt, g=GRAVITY, depth=1.0, dx=1.0, dy=1.0):
    """One midpoint (RK2) step of the linear shallow-water equations."""
    half = state + 0.5 * dt * _Derivative(state, g, depth, dx, dy)
    return state + dt * _Derivative(half, g, depth, dx, dy)


def StableTimeStep(g=GRAVITY, depth=1.0, dx=1.0, dy=1.0, cfl=0.2):
    """A CFL-guarded time step: cfl * min(dx, dy) / c with c = sqrt(g H0)."""
    wave_speed = numpy.sqrt(g * depth)
    return cfl * min(dx, dy) / wave_speed


def MakeInitialState(shape=(8, 16), seed=0, amplitude=0.1, smoothing_passes=4):
    """A smooth random height anomaly at rest: (3, H, W) with u = v = 0.

    Smoothness comes from repeated 5-point averaging of white noise -
    enough for the FD operators to resolve the field.
    """
    rng = numpy.random.default_rng(seed)
    height = rng.standard_normal(shape)
    for _ in range(smoothing_passes):
        height = (height
                  + numpy.roll(height, 1, axis=0) + numpy.roll(height, -1, axis=0)
                  + numpy.roll(height, 1, axis=1) + numpy.roll(height, -1, axis=1)) / 5.0
    height *= amplitude / max(numpy.abs(height).max(), 1e-12)
    state = numpy.zeros((3,) + tuple(shape))
    state[0] = height
    return state


def ComputeEnergy(state, g=GRAVITY, depth=1.0):
    """The quadratic invariant 0.5 * sum(g h^2 + H0 (u^2 + v^2))."""
    height, u, v = state
    return 0.5 * float(numpy.sum(g * height ** 2 + depth * (u ** 2 + v ** 2)))


def GenerateTrajectory(shape=(8, 16), steps=40, dt=None, g=GRAVITY, depth=1.0,
                       dx=1.0, dy=1.0, seed=0, amplitude=0.1, initial_state=None):
    """Integrates a smooth random initial anomaly: (T, 3, H, W) float64.

    dt defaults to StableTimeStep's CFL-guarded value; an explicit dt above
    the guard raises (the linear system would blow up).
    """
    stable = StableTimeStep(g, depth, dx, dy)
    if dt is None:
        dt = stable
    elif dt > 2.0 * stable:
        raise ValueError(
            f"dt = {dt} violates the CFL guard ({2.0 * stable:.4g} for this grid); "
            "the linear system would blow up.")
    if steps < 2:
        raise ValueError(f"steps must be >= 2 [ steps = {steps} ].")

    state = (numpy.array(initial_state, dtype=float) if initial_state is not None
             else MakeInitialState(shape, seed, amplitude))
    if state.shape != (3,) + tuple(shape):
        raise ValueError(
            f"initial_state must have shape (3, {shape[0]}, {shape[1]}); got {state.shape}.")

    trajectory = [state]
    for _ in range(steps - 1):
        state = Step(state, dt, g, depth, dx, dy)
        trajectory.append(state)
    return numpy.stack(trajectory)


def MakeStepPairs(trajectory):
    """(T, 3, H, W) trajectory -> [(state_t, state_{t+1}), ...] numpy pairs
    - the next-state training samples for a grid surrogate."""
    trajectory = numpy.asarray(trajectory, dtype=float)
    if trajectory.ndim != 4:
        raise ValueError(f"trajectory must have shape (T, C, H, W); got {trajectory.shape}.")
    return [(trajectory[t], trajectory[t + 1]) for t in range(trajectory.shape[0] - 1)]
