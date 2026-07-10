
## Convection Diffusion Application

A set of elements, conditions, strategies and solvers necessary for the solution of convection-diffusion problems.

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Readme_files/ConvectionDiffusionApplication.png" alt="Solution" style="width: 600px;"/>
</p>

The application includes tests to check the proper functioning of the application

### Features:

- A set of *Neumann* conditions:
     * Flux conditions
     * Thermal conditions

- Elements:
    * Laplacian element (both 2D/3D)
    * Eulerian convection-diffusion (both 2D/3D)
    * Convection-diffusion (both 2D/3D)
    * Convection-diffusion with change of phase (2D)
    * Explicit eulerian convection-diffusion (both 2D/3D)

- Strategies:
	* Non-linear/linear convection-diffusion strategy
	* Eulerian convection-diffusion strategy
	* Semi-Eulerian convection-diffusion strategy

- Utilities and others:
	* BFECC convection utility
	* BFECC elemental limiter convection utility
	* Convection particle
	* Face-heat utilities
	* Move particle utility
	* Pure convection tools
	* Pure convection (Crank-Nicolson) tools
