## Fluid Dynamics Application

The Fluid Dynamics Application contains the core developments in Computational Fluid Dynamics (CFD) within Kratos Multiphysics.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Kratos/assets/61457043/2dfad86e-9654-4ca2-a3e3-cb1a7baf7c63" alt="CFD_logo" width="200"/>
</p>
  
### Features:

- Stabilized FEM solvers for incompressible and compressible flow problems.

- Support for MPI parallelization (with Trilinos Application).

- Embedded formulation for flows around non-watertight or poorly defined geometries (STL) or for problems with fixed mesh and deforming domains.

- Arbitrary Lagrangian-Eulerian (ALE) formulation allows for mesh deformation during the simulation (see MeshMovingApplication).

- Support for Fluid-Structure Interaction (see FSIApplication).

- Thermally-coupled flows (with ConvectionDiffusionApplication).
