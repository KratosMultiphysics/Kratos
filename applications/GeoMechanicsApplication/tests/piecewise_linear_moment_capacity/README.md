# Piecewise Linear Moment Capacity

This test suite verifies the `PiecewiseLinearMomentCapacityConstitutiveLaw` and its integration with the
`LinearTimoshenkoBeamElement2D2N` used in the Python test cases under this folder.

## Setup

Each case models a 2 m long beam discretised with two 2-node Timoshenko beam elements. The left end is fixed firmly and the right end receives the prescribed displacement. 

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

There are four tests:

- the right end moves up linearly with time.
- the right end moves down linearly with time.
- the right end moves up then down.
- the right end moves down then up.

## Assertions

- checked displacement and bending moment at the right end.


