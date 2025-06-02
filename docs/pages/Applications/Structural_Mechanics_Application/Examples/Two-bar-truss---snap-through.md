---
title: Two bar truss   snap through
keywords: 
tags: [Two-bar-truss---snap-through.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
To demonstrate the ability of the non-linear truss element to describe geometric non-linearities the following symmetry of a two bar truss is investigated:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_truss_snap_through/symmetryTrussSystem.jpg" width="300">

_Statical System [1]_

With E = 210E09 and A = 0.01.

## Results

Two different approaches can be used to analyze the structure.

By incrementally increasing the load and solving for the residual to be zero the first limit point can be found:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_truss_snap_through/LoadCont.PNG" width="500">

_Non-linear snap through: Load-control_

Whereas both limit points can be found by incrementally increasing the displacement and solving for the residual to be zero:


<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_truss_snap_through/DispCont.PNG" width="500">

_Non-linear snap through: Displacement-control_




## References
1. Steen Krenk. Non-linear modeling and analysis of solids and structures. Cambridge
Univ. Press, 2009., pp. 28-29.