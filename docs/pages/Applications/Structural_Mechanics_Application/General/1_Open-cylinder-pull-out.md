---
title: Open cylinder pull out
keywords: 
tags: [Open-cylinder-pull-out.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
The problem considered is the static geometrically non-linear pull-out of an open cylinder with a load P = 40,000. The geometry of the cylinder is defined by: L= 10.35, R = 4.953 and h = 0.094, while the isotropic linear elastic material is characterized by E = 10.5E6 and ν = 0.3125. 

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Open_cylinder_pullout/open_cylinder_pullout_setup.png" width="600">

_Problem definition [1]_

The key displacement is the vertical deformation u_z at the point of load application as illustrated in the figure above, with the reference equilibrium path according to [1] included in the results below. 

## Results
The following deformation (exaggerated) animation of the Kratos thin quad element (mesh = 192 elements) is provided for context.

![Open cylinder pullout animation](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Open_cylinder_pullout/open_cylinder_pullout_animation.gif)

_Open cylinder pullout: Deformation of Kratos thick quad element_

The results of the test for the thin and thick triangle Kratos shell elements (mesh = 3000 elements) are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Open_cylinder_pullout/Load_displacement_curve_open_cylinder_pullout_tri.png" width="600">

_Open cylinder pullout results: triangle elements_

The results of the test for the thin and thick quadrilateral Kratos shell elements (mesh = 192 elements) are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Open_cylinder_pullout/Load_displacement_curve_open_cylinder_pullout_quad.png" width="600">

_Open cylinder pullout results: quadrilateral elements_

Both graphs above indicate all Kratos triangular and quadrilateral shell elements agree with the reference solution. 

## References
1. K.Y. Sze, X.H. Liu, and S.H. Lo. “Popular benchmark problems for geometric nonlinear analysis of shells”. In: Finite Elements in Analysis and Design 40.11 (2004), pp. 1551 –1569.