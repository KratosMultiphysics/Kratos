---
title: Shallow angle beam
keywords: 
tags: [Shallow-angle-beam.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
As the effect of bending on the shorting of the beam element is not included in the underlying co-rotational beam theory, the following example shows, how this effect can be modeled with the help of multiple beam elements. For that purpose the following shallow angled beam is investigated:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shallow_angle_beam/shallowbeamSystem.JPG" width="500">

_Statical System [1]_

With L = 10, h = 0.24, E = 210E09, ν = 0.30, A = 0.01 and Iz=Iy=IT = 3.34E-05. Whereas the load F is linearly increased.

## Results

The following deformation (exaggerated) animation of the Kratos is provided for context:

![Open cylinder pullout animation](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shallow_angle_beam/shallowAngleBeam.gif)

The results of the vertical deformation of the middle node are given in the following graph, where both the axes are scaled:

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Shallow_angle_beam/ShallowAngledBeamResult.PNG" width="500">


_Roll-up beam structure: Z-Moment_

## References
1. Steen Krenk. Non-linear modeling and analysis of solids and structures. Cambridge
Univ. Press, 2009., pp. 116-117.