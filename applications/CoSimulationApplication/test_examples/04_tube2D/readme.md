# Description of 04_tube2D

This documentation describes the test example 04_tube2D.

## Solvers

The flow solver is the pipe2D_flow solver.
The structure solver is the pipe2D_structure solver.

These solver are actual 1D, but they communicate as 2D solvers, meaning that x-, y-, z-displacemnents, pressures and x-, y-, z-traction components 
are exchanged.