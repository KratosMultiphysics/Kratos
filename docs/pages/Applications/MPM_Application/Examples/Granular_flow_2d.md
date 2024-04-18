---
title: Granular Flow 2D Validation Test
keywords: 
tags: [Granular_flow_2d.md]
sidebar: mpm_application
summary: 
---
**Author:** Bodhinanda Chandra\
**Source files:** [granular_flow_2D](https://github.com/KratosMultiphysics/Examples/tree/master/mpm/validation/granular_flow_2D/source)

## Case Specification

This is a 2D non-cohesive granular material simulation according to the experiment conducted by (Bui et al., 2008). Here, linear structured triangular elements are used to initialize the MPs and as the background mesh. The structured mesh arrangement is chosen to avoid the irregularities of the generated MP’s density, which, by further, improving the numerical solutions.

The following application of Kratos is used:
- [MPMApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication)

The problem geometry as well as the boundary conditions are sketched below:

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mpm/validation/granular_flow_2D/data/granular_flow_2D_initial.png" alt="Initial mesh and boundary conditions." width="350" />
</p>

An elasto-plastic Mohr-Coulomb plane stress constitutive law with unit thickness is considered with the following material parameters:
* Density (_&rho;_): 2650 Kg/m<sup>3</sup>
* Young's modulus (_E_):  840 kPa
* Poisson ratio (_&nu;_): 0.3
* Angle of internal friction (_&phi;_): 19.8°
* Cohesion (_c_): 0.0 kPa
* Dilatancy angle (_&psi;_): 0.0°

The time step is 0.00005 seconds, while the total simulation time is 2.0 seconds.

## Results

The problem stated above has been solved with a structured mesh with 3 material points per cell is considered with average mesh size of 0.002 m. The obtained numerical result is compared with experimental and simulation results conducted by (Bui et al., 2008) as depicted by the following figures:

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mpm/validation/granular_flow_2D/data/granular_flow_2D_results.png" alt="Obtained results and comparison." width="700" />
  
  (a.) Experiment conducted by (Bui et al., 2008), (b.) comparison of final surface configuration and failure line, (c.) simulation results of (Bui et al., 2008) by using SPH method, (d.) simulation results obtained by implicit MPM method
</p>


## References
- Bui, H. H., Fukagawa, R., Sako, K., & Ohno, S. (2008). Lagrangian meshfree particles method (SPH) for large deformation and failure flows of geomaterial using elastic-plastic soil constitutive model. International Journal for Numerical and Analytical Methods in Geomechanics, 32(12), 1537–1570. https://doi.org/10.1002/nag.688
- Chandra, B., Larese, A., Iaconeta, I., Rossi, R., Wüchner, R. (2018). Soil-Structure Interaction Simulation of Landslides Impacting a Structure Using an Implicit Material Point Method. *Accepted for publication by Proceeding of the 2nd International Conference on The Material Point Method (MPM2019)*.
