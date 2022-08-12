---
title: 2D Cavity flow problem
keywords: 
tags: [2D-Cavity-flow-problem.md]
sidebar: fluid_dynamics_application
summary: 
---


This is a 2D CFD simulation of the well known lid-driven cavity flow problem. The reference solution has been taken from  Marchi et. al. (2009), The following applications of Kratos are used:
* MeshingApplication
* ExternalSolversApplication 
* FluidDynamicsApplication

## Case Specification
A one meter square-shaped domain is considered. The fluid parameters are set such that the Reynolds number is 100, thus:
* Density (&rho;): 1.0 _Kg/m<sup>3</sup>_
* Kinematic viscosity (&nu;): 0.01 _m<sup>2</sup>/s_

The boundary conditions are:
* Bottom, left and right boundaries: No-slip
* Top boundary: _v = [1.0,0.0] m/s_
* Top left corner: _v = [0.0,0.0] m/s_
* Top right corner: _v = [1.0,0.0] m/s_
* Bottom left corner: _p = 0.0 Pa_

The time step has been set to 0.1 seconds. Regarding the total time, 10 seconds has been proved to be enough to reach a stationary solution.

## Results
The problem stated above has been solved with a structured mesh with 200x200 divisions composed by linear triangular elements. The obtained velocity field as well as a centred vertical cross section, comparing the reference solution (solved with a 1024x1024 grid) against the obtained one, are shown below. 

![Obtained velocity field.](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/2D_cavity_flow_problem/velocity_field.png)
![Cross section at _x = 0.5_ solution.](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/2D_cavity_flow_problem/vx_vertical_cut_ref2.png)

## References
Marchi, Carlos Henrique, Suero, Roberta, & Araki, Luciano Kiyoshi. (2009). The lid-driven square cavity flow: numerical solution with a 1024 x 1024 grid. Journal of the Brazilian Society of Mechanical Sciences and Engineering, 31(3), 186-198. [https://dx.doi.org/10.1590/S1678-58782009000300004](https://dx.doi.org/10.1590/S1678-58782009000300004)