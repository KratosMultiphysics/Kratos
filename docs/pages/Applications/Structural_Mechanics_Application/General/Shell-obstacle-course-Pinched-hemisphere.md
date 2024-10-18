---
title: Shell obstacle course Pinched hemisphere
keywords: 
tags: [Shell-obstacle-course-Pinched-hemisphere.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
The linear static pinched hemisphere test considers a hemispherical shell loaded with opposing point loads along it's equator. Isotropic material properties are as per the figure below. Due to symmetry only a quarter of the shell is modeled. 

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_hemisphere/pinched_hemisphere_setup.PNG" width="600">

_Problem definition [2]_

The key result is the X-displacement along one of the point loads, denoted by "u" in the following diagram, for which the reference value is "u_x" =  0.0924 [2]. 

## Results
The following X-displacement contour of the Kratos thin quad element (mesh = 285 elements) is provided for context.

![Pinched hemisphere displacement contour.](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_hemisphere/pinched_hemi_disp_contour_285elements.png)

_Pinched hemisphere results: X-displacement contour of Kratos thin quad element_

The results of the test for the thin and thick triangle Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_hemisphere/pinched_hemi_tri_results.png" width="600">

_Pinched hemisphere results: triangle elements_

The results of the test for the thin and thick quadrilateral Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shell_obstacle_course_pinched_hemisphere/pinched_hemi_quad_results.png" width="600">

_Pinched hemisphere results: quadrilateral elements_

Both graphs above indicate all Kratos triangular and quadrilateral shell elements agree with the reference solution. 

## References
1. Ted Belytschko et al. “Stress projection for membrane and shear locking in shell finite elements”. In: Computer Methods in Applied Mechanics and Engineering 51.1-3 (1985), pp. 221–258.
2. Robin Bouclier, Thomas Elguedj, and Alain Combescure. “Efficient isogeometric NURBS-based solid-shell elements: Mixed formulation and method”. In: Computer Methods in Applied Mechanics and Engineering 267 (2013), pp. 86 –110.
