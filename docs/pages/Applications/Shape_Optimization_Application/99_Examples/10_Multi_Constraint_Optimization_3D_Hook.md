---
title: Multi Constraint Optimization 3D Hook
keywords: 
tags: [10_Multi_Constraint_Optimization_3D_Hook.md]
sidebar: shape_optimization_application
summary: 
---
# Optimization of a Solid Hook
Optimization of a solid 3D Hook subjected to multiple constraints.

> **Author**: Armin Geiser
>
> **Kratos version**: 9.0

## Optimization Problem

### Objective
- Minimize mass
### Constraints
1. Strain energy of main load (LC1) &le; initial value
2. Strain energy of tip load (LC2) &le; initial value
3. No penetration of packaging (bounding) mesh

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/10_Multi_Constraint_Optimization_3D_Hook/images/hook_opt_setup.png" height="500">
</p>

## Optimization settings

- Algorithm type : `gradient_projection`
- Number of steps : `20`
- Step size : `3.0`
- Filter radius : `25.0`
- Mesh motion : `True`

## Results

### Shape Evolution
The below image shows the shape of the hook during the optimization iterations.
It can be seen how the cross section evolves towards an I-beam shape. This allows to minimize the mass is while satisfying the mechanical constraints. The shape evolution is bounded by the packaging geometry on the right side.
The color plot in the cross section shows how the internal nodes are moved by the mesh motion solver according to the shape update of the surface.

<p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/10_Multi_Constraint_Optimization_3D_Hook/images/hook_iso_mesh_color_white.gif" height="500">
</p>

### Convergence
The below plots show the evolution of objective function (i.e. mass) and the constraints over the optimization iterations.

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/10_Multi_Constraint_Optimization_3D_Hook/images/3DHookConvergencePlots.svg" height="500">
</p>

| Objective | Improvement |
| --------- | ----------- |
| mass      | `17.6%`     |

| Constraint               | Violation                  |       Remark       |
| ------------------------ | -------------------------- | :----------------: |
| strain energy: main load | `0.1%`                     |         -          |
| strain energy: tip load  | `0.05%`                    |         -          |
| packaging: bounding mesh | max nodal violation: `2.0` | `75%` of step size |

***Note:** The packaging constraint has an increasing value because more and more nodes are bounded by the constraint.*


## Source: 
[https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/10_Multi_Constraint_Optimization_3D_Hook](https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/10_Multi_Constraint_Optimization_3D_Hook)
