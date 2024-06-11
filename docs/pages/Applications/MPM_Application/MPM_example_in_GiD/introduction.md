---
title: Introduction
keywords: 
tags: [introduction.md]
sidebar: mpm_application
summary: 
---
# MPM example in GiD

## 1.Introduction
This tutorial is intended to give an overview on how to pre- and post-process a problem with the material point method (MPM) in [GiD](https://www.gidhome.com/). The general procedure how to set up a problem geometry, specify the boundary conditions, generate the mesh, calculate the solution and display the results, is explained subsequently by going through an example. For this purpose a cantilever beam under dead load is regarded. The corresponding geometric and cross-sectional properties are displayed in the picture below. For the calculation of the problem [KRATOS Multiphysics](https://github.com/KratosMultiphysics/Kratos) is called internally from GiD.

![Structural system](https://user-images.githubusercontent.com/51473791/168762544-750d2f29-6ed7-409d-8205-a6257e6a72ac.png) 
