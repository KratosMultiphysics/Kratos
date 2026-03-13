## Fluid Dynamics Application
| Left columns  | Description | Status | Authors
| ------------- | ------------| :----: | -------
| `FluidDynamicsApplication` | The Fluid Dynamics Application contains the core developments in Computational Fluid Dynamics (CFD) within Kratos Multiphysics. | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | Rubén Zorrilla (rzorrilla@cimne.upc.edu) <br /> Riccardo Rossi (rrossi@cimne.upc.edu) <br /> Jordi Cotela (jcotela@altair.com)

### 1. General features:
- Stabilized FEM solvers for incompressible, weakly-compressible and compressible flow problems.
- Support for MPI parallelization (with _MetisApplication_ and _TrilinosApplication_).
- Arbitrary Lagrangian-Eulerian (ALE) formulation allows for mesh deformation during the simulation (see _MeshMovingApplication_).
- Compatible with meshes made up with linear elements.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/use_cases/barcelona_wind/resources/Scalability.png?raw=true" alt="Wind flow over Barcelona scalability test" style="width: 600px;"/>
</p>

_Wind flow over Barcelona scalability test. More info [here](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/use_cases/barcelona_wind/README.md)._

### 2. Incompressible flows
#### Features
The simulation of viscous incompressible flows is the main capability of this application.
The application includes a variety of stabilized 2D/3D **Navier-Stokes** and **Stokes** solvers.
Limited support to 2D axisymmetric problems is also included.
Among the wide variety of stabilization techniques present in the literature, in this application the **Variational MultiScale (VMS)** (both with quasi-static and dynamic subscales), **Orthogonal SubScales (OSS)** and **Finite Increment Calculus (FIC)** methods are implemented.
All the incompressible flow elements of the application support both **Newtonian** and **non-Newtonian** (Bingham, Herschel-Bulkley) constitutive models.

A set of **boundary conditions** are included in the application. On top of the standard fixed velocity/pressure there exists the possibility to impose slip boundary conditions using MultiFreedom Constraints (MFCs) or periodic conditions using MultiPoint Constraints (MPCs).
Concerning the wall modelling, the application features linear-log and Navier-slip wall models, with the possibility to easily extend to other models.

The application also includes two different solution strategies. First one is the standard **monolithic** one in which both velocity and pressure equations are solved at once using a Newton-Raphson solver. Second one is a segregated **fractional step** strategy that accelerates the solution procedure (we note that this is only compatible with the VMS formulation).

#### Examples
- [Body-fitted 100 Re cylinder](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/validation/body_fitted_cylinder_100Re/README.md)

### 3. Weakly-compressible flows
#### Features
Similar to the described above incompressible solver, the application also includes a **VMS stabilized weakly compressible Navier-Stokes** formulation.
This solver modifies the mass conservation equation to add a slight compressibility which relates the pressure to the volume variation thanks to the inclusion of a pressure-density equation of state.
The energy equation remains uncoupled so thermal effects are assumed to be negligible.

### 4. Compressible flows
#### Features
The application includes a 2D/3D explicit compressible solver implementing a **VMS and OSS stabilized full Navier-Stokes formulations** written in **conservative variables** (momentum, density and total energy).

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

As a final note, we shall remark that at current date this solver only supports shared memory parallelism (OpenMP).

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/validation/compressible_step_woodward_colella/data/step_woodward.gif?raw=true" alt="Woodward and Colella's Mach 3 step density field." style="width: 600px;"/>
</p>

_Woodward and Colella's Mach 3 step density field._

#### Examples
- [Transonic flow around a NACA0012 profile](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_naca_0012_Ma_0.8/README.md)
- [Multistage transonic flow around a NACA0012 profile](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/multistage_compressible_naca_0012_Ma_0.8/README.md)
- [Transonic flow around a NACA0012 profile at a 3&deg; angle](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_naca_0012_Ma_0.8_aoa_3/README.md)
- [Sod shock tube](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_sod_shock_tube/README.md)
- [Supersonic flow in Woodward and Colella's Mach 3 step](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_step_woodward_colella/README.md)
- [Supersonic flow over a wedge](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/compressible_Wedge/README.md)

### 5. Unfitted mesh methods
#### Features
The embedded solver allows the resolution of problems with **unffitted boundaries**, including flows around volumetric and volumeless (i.e. shell-like) bodies.
Starting from a distance field, either analytical or obtained with any of the levelset algorithms in _KratosCore_, the embedded solver uses a **Cut-FEM** approach to solve the problem.
This approach only supports simplicial meshes. (linear triangle and tetrahedron).
This solver, which can be used in combination with all the formulations described in the incompressible flow section, makes possible to efficiently solve flows around non-watertight or poorly defined geometries (e.g. STL) as well as cases involving arbitrary large boundary displacements and rotations.

Current research on this topic include the development of **Shifted Boundary Method (SBM)** solvers.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/validation/embedded_moving_cylinder/data/embedded_moving_cylinder_p.gif?raw=true" alt="Embedded moving cylinder velocity field [m/s]." style="width: 600px;"/>
</p>

_Embedded moving cylinder example velocity field._

#### Examples
- [Embedded moving cylinder](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/validation/embedded_moving_cylinder/README.md)

### 6. Two-phase flows
#### Features
The _FluidDynamicsApplication_ includes a solver for the resolution of **biphasic (Newtonian-air) viscous incompressible flows**.
This solver uses a **levelset based approach** which combines the **implicit** fluid solver with a convection and redistancing ones (see the _KratosCore_ for more information about these).
The solver is able to account for the pressure discontinuities thanks to a local enrichment plus an element-by-element static condensation, which avoids the need to reform the sparse matrix graph at each time step.
Besides, the solver is also equipped with a strategy to revert the mass losses introduced by the levelset approach.

#### Examples
- [Two-fluids dam break scenario](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/validation/two_fluid_dam_break/README.md)
- [Two-fluids wave propagation](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/validation/two_fluid_wave/README.md)

### 7. Multiscale modelling
The application also includes limited support for the multiscale modelling following the **Representative Volume Element (RVE)** approach.

### 8. Multiphysics problems
#### Features
The _FluidDynamicsApplication_ can be coupled with other applications to solve multiphysics problems such as **Fluid-Structure Interaction (FSI)** (see _FSIApplication_) or **thermally-coupled** flows (buoyancy and Conjugate Heat Transfer (CHT)) (see _ConvectionDiffusionApplication_).

#### Examples
Conjugate Heat Transfer:
- [Cylinder cooling Re = 100 and Pr = 2](https://github.com/KratosMultiphysics/Examples/blob/master/conjugate_heat_transfer/validation/cylinder_cooling_Re100_Pr2/README.md)

Fluid-Structure Interaction:
- [FSI lid driven cavity](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_structure_interaction/validation/fsi_lid_driven_cavity/README.md)
- [Mixer with flexible blades (embedded)](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_structure_interaction/validation/embedded_fsi_mixer_Y/README.md)
- [Mok benchmark](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_structure_interaction/validation/fsi_mok/README.md)
- [Mok benchmark (embedded)](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_structure_interaction/validation/embedded_fsi_mok/README.md)
- [Turek benchmark - FSI2](https://github.com/KratosMultiphysics/Examples/blob/master/fluid_structure_interaction/validation/fsi_turek_FSI2/README.md)
