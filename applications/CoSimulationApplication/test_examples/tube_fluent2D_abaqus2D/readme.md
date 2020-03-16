# Description of tube2D_fluent_abaqus

This documentation describes the test example tube2D_fluent_abaqus.

## Solvers

The flow solver is fluent, used to solve an axisymmetric representation of a tube. The axial direction is along the X-axis,
the radial direction along the Y-axis.

The structure solver is abaqus, used to solve an axisymmetric representation of a tube. The axial direction is along the Y-axis,
the radial direction along the X-axis.

A permutation mapper is introduced between the solvers, to match the coordinate frames.
