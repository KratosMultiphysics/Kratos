---
title: Shape Update Optimization Stanford Bunny
keywords: 
tags: [Shape_Update_Optimization_Stanford_Bunny.md]
sidebar: 
summary: 
---
# Smooth surface wrapping

A pure geometric optimization problem of a wrapping surface smoothly around a complex object, which in this example is a Stanford bunny.

> **Author**: Armin Geiser
>
> **Kratos version**: 9.0

## Optimization Problem

### Objective
- Shape update maximization

### Constraints
- No penetration of packaging - bounding mesh (*Standford bunny*)

  <p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/11_Shape_Update_Optimization_Stanford_Bunny/images/bunny_opt_setup.png" height="500">
  </p>

## Optimization settings
- Algorithm type : `gradient_projection`
- Number of steps : `150`
- Step size : `0.001`
- Filter radius : `0.015`
- Mesh motion : `False`

## Results

### Shape Evolution
The below image shows the shape evolution of the wrapping surface during the optimization iterations. The small sphere *grows* inside the bounding geometry. Using the relatively large filter radius, a smooth local optimum is found. Details smaller than the filter radius (e.g. the ears) are not captured.

<p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/11_Shape_Update_Optimization_Stanford_Bunny/images/bunny_results_smallSphere.gif" height="400">
</p>


## Alternative problem setup

The bounding geometry of the Stanford bunny can be approached (wrapped) from the opposite site as well. For such a *shrinking* optimization we start from a large sphere. Compared to the example described above, the starting geometry is exchanged, the feasible side for the bounding geometry is switched and the objective is swtiched as well - see the setup in the [shrink folder](shrink) for details.

  <p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/11_Shape_Update_Optimization_Stanford_Bunny/images/bunny_results_largeSphere.gif" height="500">
  </p>



## Source: 
[https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/11_Shape_Update_Optimization_Stanford_Bunny](https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/11_Shape_Update_Optimization_Stanford_Bunny)
