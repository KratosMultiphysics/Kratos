 
## Structural Mechanics Application  
  
The Structural Mechanics Application contains a series of structural elements, as well as solid elements, constitutive laws and the corresponding strategies and solvers within Kratos Multiphysics.  
 
<p align="center"> 
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/beam_roll_up/data/rollup.gif" alt="Solution" style="width: 600px;"/> 
</p> 
  
The application includes tests to check the proper functioning of the application 
  
### Features:  
  
- A set of *Neumann* conditions:
     * Point loads (loads applied directly on the nodes)
     * Point moment (a discret moment applied directly on the nodes) 
     * Line load (a distributed load applied over a line)
     * Surface load (a distributed load applied over a face)
     * A simple point contact conditions based on the distance
     
- Solid elements:
    * Small displacement elements
    * Total Lagrangian elements
    * Updated Lagrangian elements
    
- Structural elements:
    * Zero-dimensional elements :
        * Nodal concentrated element (both 2D/3D). Includes nodal damping, nodal mass and nodal stiffness
    * Uni-dimensional elements :
        * Truss element (3D)
       	* Corrotational beam element (both 2D/3D)
    * Two-dimensional elements :
        * Membrane (regular and pre-stress)
        * Isotropic shell element
        * Thin shell (Quadrilateral and triangular)
       	* Thick shell (Quadrilateral and triangular)
       	
- Constitutive laws: 
	* Isotropic laws (Plane strain, plane stress and 3D)
	* Orthotropic law (Plane stress)
	* Hyperelastic laws:
		* Neo-Hookean
		* Kirchhoff
		
- Strategies:
	* Formfinding strategies
	* Eigensolver strategy
	* Harmonic analysis strategies
	* Arc-length strategy
	
- Schemes:
	* Relaxation scheme
	* Eigen solver scheme
	
- Convergence criteria:
	* For displacement and other DoF
	* For displacement and rotation
	
- Builder and solver:    
	* Multi-point Constraint builder and solver
	
- Utilities and others:
	* A process to post-process eigenvalues
	* A GiDIO utility for eigen values
