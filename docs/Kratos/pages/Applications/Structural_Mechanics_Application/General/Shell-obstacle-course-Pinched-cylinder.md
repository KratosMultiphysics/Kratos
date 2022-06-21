---
title: Shell obstacle course Pinched cylinder
keywords: 
tags: [Shell-obstacle-course-Pinched-cylinder.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
The linear static pinched cylinder test considers a cylindrical shell fixed by rigid diaphragms at it's axial ends. The loading consists of two opposing compressive point loads at the center of the shell. Isotropic material properties are as per the figure below. Due to symmetry only an eighth of the shell is modeled. 

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_cylinder/pinched_cylinder_setup.PNG" width="600">

_Problem definition [2]_

The key result is the vertical displacement under the point load, denoted by "u" in the diagram above, for which the reference value is u_z =  1.8248E-5 [2]. 

## Results
The following Z-displacement contour of the Kratos thin quad element (mesh = 384 elements) is provided for context.

![Pinched cylinder displacement contour.](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_cylinder/pinched_cylinder_disp_contour_384elements.png)

_Pinched cylinder results: Z-displacement contour of Kratos thin quad element_

The results of the test for the thin and thick triangle Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_cylinder/pinched_cyl_structured_tri_results.png" width="600">

_Pinched cylinder results: triangle elements_

The results of the test for the thin and thick quadrilateral Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_cylinder/pinched_cyl_structured_quad_results.png" width="600">

_Pinched cylinder results: quadrilateral elements_

Both graphs above indicate all Kratos triangular and quadrilateral shell elements agree with the reference solution. 

## References
1. Ted Belytschko et al. “Stress projection for membrane and shear locking in shell finite elements”. In: Computer Methods in Applied Mechanics and Engineering 51.1-3 (1985), pp. 221–258.
2. Robin Bouclier, Thomas Elguedj, and Alain Combescure. “Efficient isogeometric NURBS-based solid-shell elements: Mixed formulation and method”. In: Computer Methods in Applied Mechanics and Engineering 267 (2013), pp. 86 –110.
