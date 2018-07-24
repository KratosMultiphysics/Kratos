## Fluid Dynamics Application

The Fluid Dynamics Application contains the core developments in Computational Fluid Dynamics (CFD) within Kratos Multiphysics.

### Features:

- Stabilized FEM solvers for incompressible flow problems.

- Support for MPI parallelization (with Trilinos Application).

- Embedded formulation for flows around non-watertight or poorly defined geometries (STL) or for problems with fixed mesh and deforming domains.

- Arbitrary Lagrangian-Eulerian (ALE) formulation allows for mesh deformation during the simulation (see MeshMovingApplication).

- Support for Fluid-Structure Interaction (see FSIapplication).

- Thermally-coupled flows (with ConvectionDiffusionApplication).
