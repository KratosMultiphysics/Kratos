---
title: Cylinder on inclined plane 2D - comparison between analytical and numerical solution with MPM
keywords: 
tags: [Cylinder_on_inclined_plane.md]
sidebar: mpm_application
summary: 
---
**Author:** Philip Franz\
**Source files:** [cylinder_on_inclined_plane_2D](https://github.com/KratosMultiphysics/Examples/tree/master/mpm/validation/cylinder_on_inclined_plane/source)

## Case Specification

This is a 2D simulation of a cylinder on an inclined plane. A rotating as well as a frictionless sliding behaviour of the cylinder are regarded subsequently. The simulation is set up according to section 4.5.2 of (Iaconeta, 2019). 
Linear, unstructured, triangular elements with a size of 0.01m are used to initialize the MPs. Three MPs per cell are considered. For the backgroundmesh linear, unstructured, triangular elements with a size of 0.02m are used.
However, in contrast to section 4.5.2 of (Iaconeta, 2019), the inclined plane is modelled by a line with unstructured elements with size 0.01m. On that line a non conforming Dirichlet boundary condition is imposed by using the penalty method based on (Chandra et al., 2021).  

The following applications of Kratos are used:
- [MPMApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication)
- [LinearSolversApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/LinearSolversApplication)

The problem geometry as well as the boundary conditions are sketched below. The non conforming boundary condition is respresented by the copper coloured line.

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mpm/validation/cylinder_on_inclined_plane/data/Cylinder_on_inclined_plane_with_grid_detail.png" alt="Initial geometry and boundary conditions." alt="Initial geometry and boundary conditions."  width="1400" />
</p>

A hyper elastic Neo Hookean Plane strain (2D) constitutive law with unit thickness is considered with the following material parameters:
* Density (_&rho;_): 7800 Kg/m<sup>3</sup>
* Young's modulus (_E_):  200 MPa
* Poisson ratio (_&nu;_): 0.3

The time step is 0.001 seconds; the total simulation time is 1.0 seconds. The angle (_&alpha;_) of the inclined plane is 60°. The penalty-factor is 1e13. 

The contact between cylinder and inclined plane is modelled with the option "contact" (see line 53, file *ProjectParameters_contact.json*) in the first and with "slip" in the second case, based on (Chandra et al., 2021). Choosing "contact" leads to a rolling behaviour of the cylinder; "slip" to frictionless sliding.

## Results
The analytical and numerical solution for the displacement function of the respective case of the above stated problem are compared afterwards:

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mpm/validation/cylinder_on_inclined_plane/data/comparison_analytical_numerical_disp.png" alt="Initial geometry and boundary conditions." width="1600" />
</p>


The left image displays the rolling cylinder - modelled with option "contact". The right one shows the sliding cylinder (frictionless) - modelled with option "slip".

<p align="center">
  <img alt="Light" src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mpm/validation/cylinder_on_inclined_plane/data/rolling cylinder gif.gif" width="45%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img alt="Dark" src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mpm/validation/cylinder_on_inclined_plane/data/sliding cylinder gif.gif" width="46%">
</p>

 
## References
- Iaconeta, I. (2019). *Discrete-continuum hybrid modelling of flowing and static regimes.* (Ph.D. thesis). Universitat politècnica de Catalunya - Barcelona tech 
- Chandra, B., Singer, V., Teschemacher, T., Wüchner, R., Larese, A. (2021) *Nonconforming Dirichlet boundary conditions in implicit material point method by means of penalty augmentation*. Acta Geotech. 16, 2315–2335. https://doi.org/10.1007/s11440-020-01123-3 
