---
title: Roll up beam cantilever
keywords: 
tags: [Roll-up-beam-cantilever.md]
sidebar: structural_mechanics_application
summary: 
---

## Problem definition
To show the ability of the co-rotational beam element to handle large rotations the following example describes a roll up beam structure. For that purpose a cantilever of the length L = 1.00 is loaded with a linearly increased z-moment at its tip. The structure is discretized with 100 elements, E = 210E09, ν = 0.30, A = 0.01 and Iz=Iy=IT = 0.00001.

<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Roll_up_beam_cantilever/RollUpSystem.JPG" width="500">

_Statical System [1]_

## Results

The beam will form one total circal as soon as the curvature is κ = 1/R = 2π/L = 6.283 and therefore the z-moment equals M_z = κ*EI_z = 1.32E07.

![Open cylinder pullout animation](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Roll_up_beam_cantilever/rollup.gif)

_Roll-up beam struture: Z-Moment_

## References
1. Steen Krenk. Non-linear modeling and analysis of solids and structures. Cambridge
Univ. Press, 2009., pp. 114-115.