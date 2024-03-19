---
title: Non linear cantilever beam
keywords: 
tags: [Non-linear-cantilever-beam.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
To demonstrate the ability of the co-rotational beam element to handle very large deformation a cantilever with L = 2.00 and a linearly increased transverse load is investigated:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_cantilever_beam/beamnonlinearcantileverSys.jpg" width="400">

_Statical System [1]_

With E = 210E09, ν = 0.30, A = 0.01 and Iz=Iy=IT = 0.00001.

## Results

By discretizing the total cantilever with 20 elements the following results, with respect to the three degrees of freedom, are obtained:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_cantilever_beam/CantV.PNG" width="500">

_Non-linear cantilever: Vertical displacement (20 elements)_

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_cantilever_beam/CantU.PNG" width="500">

_Non-linear cantilever: Horizontal displacement (20 elements)_

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_cantilever_beam/Cantphi.PNG" width="500">

_Non-linear cantilever: z-rotation (20 elements)_


A convergence study is done with respect to the horizontal dispalcement:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Non_linear_cantilever_beam/uConverg.PNG" width="500">

_Non-linear cantilever: Convergence study - horizontal displacement_



## References
1. Steen Krenk. Non-linear modeling and analysis of solids and structures. Cambridge
Univ. Press, 2009., pp. 115-116.