# Piecewise Linear Moment Capacity

This test suite verifies the `PiecewiseLinearMomentCapacityConstitutiveLaw` and its integration with the
`LinearTimoshenkoBeamElement2D2N` used in the Python test cases under this folder.

## Setup

Each case models a 2 m long beam discretised with two 2-node Timoshenko beam elements. The left end is fixed and the right end receives the prescribed displacement/load table. 

## Material properties

The tests use the property set defined in `common/MaterialParameters.json`. Current values in that file are:

- `YOUNG_MODULUS`: 200.0
- `POISSON_RATIO`: 0.3
- `UNRELOAD_MODULUS`: 200.0
- `THICKNESS`: 0.1
- `THICKNESS_EFFECTIVE_Y`: 0.1
- `KAPPA_PIECEWISE_LINEAR_LAW`: [1.0e-02, 0.5, 1.0]
- `MOMENTUM_PIECEWISE_LINEAR_LAW`: [0.1, 1.0, 1.5]

## Tests (cases)

The folder contains four case subfolders. Each subfolder includes the MDPA and (generated) GiD output files:

- `tension/` — prescribed displacement sequence defined in the MDPA table (file: `tension/piecewise_linear_moment_capacity.mdpa`).
- `compression/` — compression-driven case (file: `compression/piecewise_linear_moment_capacity.mdpa`).
- `tension_compression/` — a load history that combines tension then compression.
- `compression_tension/` — the reversed load history.

## Assertions

	- Element `1` gauss-point `AXIAL_FORCE` from the `.post.res` file and compare it to expected arrays.
	- Node `2` `DISPLACEMENT` x-component (mid node) and compare to expected arrays.


