---
title: Shell obstacle course Scordelis Lo roof
keywords: 
tags: [Shell-obstacle-course-Scordelis-Lo-roof.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
The linear static Scordelis-Lo roof test considers a section of cylindrical shell fixed by rigid diaphragms at it's axial ends subject to a pseudo-gravity distributed load of a magnitude 90. Isotropic material properties are as per the figure below. Due to symmetry, only a quarter of the shell is modeled. 

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_scordelis_lo_roof/scordelis_lo_problem_setup.PNG" width="600">

_Problem definition [2]_

The key result is the vertical displacement of the lateral side at the midpoint, denoted by "u_z" in the diagram above, for which the reference value is "u_z" = 0.3024 [2].

## Results
The following Z-displacement contour of the Kratos thin quad element (mesh = 192 elements) is provided for context.

![Scordelis Lo displacement contour.](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_scordelis_lo_roof/scordelis_lo_contour.png)

_Scordelis-Lo roof results: Z-displacement contour of Kratos thin quad element_

The results of the test for the thin and thick triangle Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_scordelis_lo_roof/scordelis_structured_tri_results_kratos.png" width="600">

_Scordelis-Lo roof results: triangle elements_

The results of the test for the thin and thick quadrilateral Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_scordelis_lo_roof/scordelis_structured_quad_results_kratos.png" width="600">

_Scordelis-Lo roof results: quadrilateral elements_

Both graphs above indicate all Kratos triangular and quadrilateral shell elements agree with the reference solution. 

## References
1. Ted Belytschko et al. “Stress projection for membrane and shear locking in shell finite elements”. In: Computer Methods in Applied Mechanics and Engineering 51.1-3 (1985), pp. 221–258.
2. Robin Bouclier, Thomas Elguedj, and Alain Combescure. “Efficient isogeometric NURBS-based solid-shell elements: Mixed formulation and method”. In: Computer Methods in Applied Mechanics and Engineering 267 (2013), pp. 86 –110.
