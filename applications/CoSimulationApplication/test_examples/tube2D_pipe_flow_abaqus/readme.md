# Description of pipe_flow_abaqus

This documentation describes the test example pipe_flow_abaqus.

## Solvers

The flow solver is pipe_flow, a 1D solver which communicates as a 2D solver. The axial direction is along the Z-axis,
the radial direction along the Y-axis.

The structure solver is abaqus, which solves an axisymmetric case. The axial direction is along the Y-axis,
the radial direction along the X-axis.

A permutation mapper is introduced between the solvers, to match the coordinate frames.
