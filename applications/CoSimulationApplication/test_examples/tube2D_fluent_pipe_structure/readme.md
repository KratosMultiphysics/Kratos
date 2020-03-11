# Description of fluent_pipe_structure

This documentation describes the test example fluent_pipe_structure.

## Solvers

The flow solver is fluent, used to solve an axisymmetric representation of a tube. The axial direction is along the X-axis,
the radial direction along the Y-axis.

The structure solver is the 2D pipe_structure_inert solver. The axial direction is along the Z-axis,
the radial direction along the Y-axis.

A permutation mapper is introduced between the solvers, to match the coordinate frames.
