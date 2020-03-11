# Description of tube3D_fluent_abaqus

This documentation describes the test example tube3D_fluent_abaqus.

## Solvers

The flow solver is fluent. The axial direction is along the X-axis,
the radial direction along the Y-axis.

The structure solver is abaqus. The axial direction is along the Y-axis,
the radial direction along the X-axis.

A permutation mapper is introduced between the solvers, to match the coordinate frames.
