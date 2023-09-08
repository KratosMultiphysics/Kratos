## Fluid Dynamics Application

The Fluid Dynamics Application contains the core developments in Computational Fluid Dynamics (CFD) within Kratos Multiphysics.

### General features:

- Stabilized FEM solvers for incompressible and compressible flow problems.

- Support for MPI parallelization (with Trilinos Application).

- Embedded formulation for flows around non-watertight or poorly defined geometries (STL) or for problems with fixed mesh and deforming domains.

- Arbitrary Lagrangian-Eulerian (ALE) formulation allows for mesh deformation during the simulation (see MeshMovingApplication).

- Support for Fluid-Structure Interaction (see FSIApplication).

- Thermally-coupled flows (with ConvectionDiffusionApplication).

### Incompressible flow



### Weakly-compressible flow
Similar to the described above incompressible solver, the application also includes a **VMS stabilized weakly compressible Navier-Stokes** formulation.
This solver modifies the mass conservation equation to add a slight compressibility which relates the pressure to the volume variation thanks to the inclusion of a pressure-density equation of state.
The energy equation remains uncoupled so thermal effects are assummed to be negligible.

### Compressible flow
#### Features
The application includes a 2D/3D explicit compressible solver implementing a **VMS stabilized full Navier-Stokes formulation** written in **conservative variables** (momentum, density and total energy).

A set of **explicit strategies** can be used
- Forward Euler
- Midpoint rule
- 3rd order Total Variational Diminishing Runge-Kutta (RK3-TVD)
- 4th order Runge-Kutta (RK4)

Two different **shock capturing** techniques are provided
- Physics-based shock capturing
- Entropy-based shock capturing

This solver can be combined in a multistage fashion with the ones in the [_CompressiblePotentialFlowApplication_](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/CompressiblePotentialFlowApplication).
By doing so, the potential solution can be used as initial condition to ease and accelerate the convergence of the full Navier-Stokes simulation.

#### Examples
- [Transonic flow around a NACA0012 profile](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_naca_0012_Ma_0.8/README.md)
- [Multistage transonic flow around a NACA0012 profile](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/multistage_compressible_naca_0012_Ma_0.8/README.md)
- [Transonic flow around a NACA0012 profile at a 3&deg; angle](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_naca_0012_Ma_0.8_aoa_3/README.md)
- [Sod shock tube](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_sod_shock_tube/README.md)
- [Supersonic flow in Woodward and Colella's Mach 3 step](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_step_woodward_colella/README.md)
- [Supersonic flow over a wedge](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_Wedge/README.md)

### Two-phase flow

### Unfitted mesh methods


### Multiscale modelling


### Coupling
- Buoyancy
- CHT
- FSI


