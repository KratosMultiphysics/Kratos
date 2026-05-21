---
title: Composite hinged cylindrical roof snapthrough
keywords: 
tags: [Composite-hinged-cylindrical-roof-snapthrough.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
In the static geometrically non-linear composite hinged cylindrical roof snapthrough test, the [isotropic hinged cylindrical roof snapthrough](Hinged-cylindrical-roof-snapthrough)
test is identically repeated with an orthotropic laminate material. The composite laminate of total thickness h = 12.7mm is composed of orthotropic laminae with the following properties:

* E_1 = 3.3 GPa
* E_2 = 1.1 GPa
* G_12 = G_13 = G_23 = 0.66 GPa
* ν_12 = 0.25

The material is arranged in a [90/0/90] stacking sequence [1].

The key result is the vertical point displacement under the point load P in the diagram above, for which the reference equilibrium path according to [1] is plot in the results below. 

## Results
The following deformation (exaggerated) animation of the Kratos thin quad element (mesh = 256 elements) is provided for context.

![Composite hinged cylindrical roof snapthrough animation](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Composite_hinged_cylindrical_roof_snapthrough/hinged_cylindrical_roof_animation.gif)

_Composite hinged cylindrical roof snapthrough: Deformation of Kratos thin quad element_

The results of the test for the thin quadrilateral and thick triangle Kratos shell elements are presented below.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Composite_hinged_cylindrical_roof_snapthrough/Load_displacement_curve_composite_hinged_cylindrical_roof.png" width="600">

_Composite hinged cylindrical roof snapthrough results: thin quadrilateral and thick triangle elements_

The graph above indicate the thick triangular and thin quadrilateral Kratos shell elements agree with the reference solution. 

## References
1. K.Y. Sze, X.H. Liu, and S.H. Lo. “Popular benchmark problems for geometric nonlinear analysis of shells”. In: Finite Elements in Analysis and Design 40.11 (2004), pp. 1551 –1569.