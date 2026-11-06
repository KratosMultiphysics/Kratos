"""Synthetic Lennard-Jones molecular-dynamics reference data for particle
surrogates.

A pure-numpy velocity-Verlet integrator of N identical atoms in a periodic
box under the truncated-and-shifted Lennard-Jones potential, in reduced
units (sigma = epsilon = mass = 1):

    V(r) = 4 [ (1/r)^12 - (1/r)^6 ] - V(r_c)      for r < r_c, else 0

with the minimum-image convention for every pair. It exists to feed the
Lennard-Jones GNN recipe - NVIDIA's own examples/molecular_dynamics/
lennard_jones trains the generic MeshGraphNet on (positions -> forces)
frames from OpenMM - when neither OpenMM nor a compiled DEMApplication is
at hand, exactly as shallow_water_reference stands in for
ShallowWaterApplication.

Positions are returned UNWRAPPED (continuous trajectories): the
minimum-image convention makes the forces periodic regardless, and
CreateParticleTrajectoryDataset's finite differences never see a boundary
jump. Because velocity Verlet's position update is the Stoermer-Verlet
recurrence x_{t+1} = 2 x_t - x_{t-1} + dt^2 a_t, the dataset's
central-difference acceleration targets equal the forces at x_t EXACTLY
(mass 1), which is what makes the recipe checkable against this module.

numpy-only and eagerly importable: no torch, no physicsnemo.
"""

import numpy

DEFAULT_CUTOFF = 2.5


def _ResolveCutoff(cutoff, box):
    """The cutoff, defaulting to the usual 2.5 sigma or half the smallest
    box length, whichever is smaller (the minimum image is unambiguous
    only up to half the box)."""
    limit = 0.5 * float(box.min())
    if cutoff is None:
        return min(DEFAULT_CUTOFF, limit)
    cutoff = float(cutoff)
    if cutoff <= 0.0 or cutoff > limit:
        raise ValueError(
            f"cutoff must lie in (0, box/2] = (0, {limit}] for the minimum image to be "
            f"unambiguous; got {cutoff}.")
    return cutoff


def _Box(box_size):
    box = numpy.asarray(box_size, dtype=numpy.float64).reshape(-1)
    if box.size == 1:
        box = numpy.repeat(box, 3)
    if box.size != 3 or numpy.any(box <= 0.0):
        raise ValueError(
            f"box_size must be one positive length or three [Lx, Ly, Lz]; got {box_size}.")
    return box


def MinimumImage(displacements, box_size):
    """The shortest periodic image of (..., 3) displacements."""
    box = _Box(box_size)
    displacements = numpy.asarray(displacements, dtype=numpy.float64)
    return displacements - box * numpy.round(displacements / box)


def ComputeForcesAndPotential(positions, box_size, cutoff=None):
    """Per-atom forces (N, 3) and per-atom potential energies (N,).

    Exact O(N^2) minimum-image evaluation - the reference, not a fast
    path. The potential is truncated and shifted at the cutoff (None: 2.5
    sigma, or half the box if that is smaller), so it is continuous and
    the total energy is a conserved quantity of the integrator; each
    atom's energy is half its pair sum, so the per-atom values sum to the
    total potential energy.
    """
    box = _Box(box_size)
    positions = numpy.asarray(positions, dtype=numpy.float64).reshape(-1, 3)
    cutoff = _ResolveCutoff(cutoff, box)

    deltas = MinimumImage(positions[:, None, :] - positions[None, :, :], box)  # r_i - r_j
    r2 = numpy.einsum("ijk,ijk->ij", deltas, deltas)
    numpy.fill_diagonal(r2, numpy.inf)
    inside = r2 < cutoff ** 2
    inv_r2 = numpy.where(inside, 1.0 / numpy.where(inside, r2, 1.0), 0.0)
    inv_r6 = inv_r2 ** 3

    pair_force = 24.0 * (2.0 * inv_r6 ** 2 - inv_r6) * inv_r2  # F_ij / |r_ij| along r_i - r_j
    forces = numpy.einsum("ij,ijk->ik", pair_force, deltas)
    shift = 4.0 * (cutoff ** -12 - cutoff ** -6)
    pair_energy = numpy.where(inside, 4.0 * (inv_r6 ** 2 - inv_r6) - shift, 0.0)
    potential = 0.5 * pair_energy.sum(axis=1)
    return forces, potential


def Step(positions, velocities, dt, box_size, cutoff=None, forces=None):
    """One velocity-Verlet step.

    Returns (positions, velocities, forces) with the forces evaluated at
    the NEW positions, so a caller can hand them back in and pay one force
    evaluation per step.
    """
    if dt <= 0.0:
        raise ValueError(f"dt must be > 0 [ dt = {dt} ].")
    positions = numpy.asarray(positions, dtype=numpy.float64)
    velocities = numpy.asarray(velocities, dtype=numpy.float64)
    if forces is None:
        forces, _ = ComputeForcesAndPotential(positions, box_size, cutoff)
    half = velocities + 0.5 * dt * forces
    positions = positions + dt * half
    forces, _ = ComputeForcesAndPotential(positions, box_size, cutoff)
    velocities = half + 0.5 * dt * forces
    return positions, velocities, forces


def MakeLattice(atoms_per_side, box_size):
    """A simple cubic lattice of atoms_per_side^3 atoms filling the box."""
    box = _Box(box_size)
    n = int(atoms_per_side)
    if n < 1:
        raise ValueError(f"atoms_per_side must be >= 1 [ atoms_per_side = {atoms_per_side} ].")
    axis = numpy.arange(n) / n
    grid = numpy.stack(numpy.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    return (grid + 0.5 / n) * box


def MakeInitialState(atoms_per_side=4, box_size=None, temperature=0.5, seed=0):
    """Lattice positions plus Maxwell-Boltzmann velocities with zero total
    momentum: (positions (N, 3), velocities (N, 3), box (3,)).

    The default box gives a lattice spacing of 1.5 sigma (a dilute fluid
    whose nearest neighbours sit on the attractive side of the well), so
    the forces are neither zero nor stiff at the default time step.
    """
    box = _Box(1.5 * atoms_per_side if box_size is None else box_size)
    positions = MakeLattice(atoms_per_side, box)
    rng = numpy.random.default_rng(seed)
    velocities = rng.standard_normal(positions.shape) * numpy.sqrt(max(float(temperature), 0.0))
    velocities -= velocities.mean(axis=0)
    return positions, velocities, box


def ComputeEnergy(positions, velocities, box_size, cutoff=None):
    """Total energy: kinetic plus the shifted pair potential."""
    velocities = numpy.asarray(velocities, dtype=numpy.float64)
    _, potential = ComputeForcesAndPotential(positions, box_size, cutoff)
    return 0.5 * float(numpy.sum(velocities ** 2)) + float(potential.sum())


def GenerateTrajectory(atoms_per_side=4, steps=50, dt=0.005, box_size=None,
                       cutoff=None, temperature=0.5, seed=0):
    """Integrates a lattice-plus-thermal-velocities initial state.

    Returns a dict:
        "positions"  (T, N, 3) UNWRAPPED float64,
        "velocities" (T, N, 3),
        "forces"     (T, N, 3)  at each frame's positions,
        "potential"  (T, N)     per-atom potential energies,
        "box_size"   (3,),
        "dt", "cutoff" (the resolved value).
    """
    if steps < 2:
        raise ValueError(f"steps must be >= 2 [ steps = {steps} ].")
    positions, velocities, box = MakeInitialState(atoms_per_side, box_size, temperature, seed)
    cutoff = _ResolveCutoff(cutoff, box)
    forces, potential = ComputeForcesAndPotential(positions, box, cutoff)
    frames = [(positions, velocities, forces, potential)]
    for _ in range(steps - 1):
        positions, velocities, forces = Step(positions, velocities, dt, box, cutoff, forces)
        _, potential = ComputeForcesAndPotential(positions, box, cutoff)
        frames.append((positions, velocities, forces, potential))
    return {
        "positions": numpy.stack([f[0] for f in frames]),
        "velocities": numpy.stack([f[1] for f in frames]),
        "forces": numpy.stack([f[2] for f in frames]),
        "potential": numpy.stack([f[3] for f in frames]),
        "box_size": box,
        "dt": float(dt),
        "cutoff": float(cutoff),
    }
