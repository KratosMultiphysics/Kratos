 
## Convection Diffusion Application  
  
The Convection DIffusion Application contains a series of elements and conditions and the corresponding strategies and solvers within Kratos Multiphysics necesaries in order to simulate a convection-diffusion problem.  
 
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
