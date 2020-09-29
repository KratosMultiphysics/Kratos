## Geo-Mechanics Application

The Geo-Mechanics Application contains features needed for common geotechnical/geomechanical applications within Kratos Multiphysics.

### Features:

- UPw small displacement element for saturated porous media (with
equal order interpolation, unstable under incompressible-undrained
conditions)

- Stable UPw small displacement element for saturated porous media
(with higher order interpolation for displacements)

- FIC-Stabilized UPw small displacement element for saturated porous media
(with equal order interpolation for displacements)

- UPw Quasi-zero-thickness interface elements for defining cracks and
joints


### How to use MPI in Geo-Mechanics Application

Make sure that the following lines are properly set in the configure.sh file:

-DEXTERNAL_SOLVERS_APPLICATION=ON        \
-DSTRUCTURAL_MECHANICS_APPLICATION=ON    \
-DGEO_MECHANICS_APPLICATION=ON           \
-DINSTALL_EMBEDDED_PYTHON=ON             \



*Note*: For the moment, MPI does not work.
