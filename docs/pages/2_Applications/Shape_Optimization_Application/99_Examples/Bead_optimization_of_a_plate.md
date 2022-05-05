---
title: Bead optimization of a plate
keywords: 
tags: [Bead_optimization_of_a_plate.md]
sidebar: shape_optimization_application
summary: 
---
# Bead Optimization of a Plate

A bead optimization problem to minimize the strain energy of a square plate - supported at the corners and subject to a single force in the center.

> **Author**: Armin Geiser
>
> **Kratos version**: 9.0

## Optimization Problem

### Objective
- Minimize strain energy

  <p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_SetupwithBC.png">
  </p>

## Optimization settings
- Algorithm type : `bead_optimization`
- Number of steps : `300`
- Step size : `0.25`
- Filter radius : `0.075`
- Mesh motion : `False`

## Results

### Shape Evolution
The below images shows the shape evolution of the bead pattern on a plate during the optimization iterations. 

<p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_result.gif" width="800">
</p>

### Comparison with manual solutions

In the table below, a comparison with manual solutions shows that the bead pattern obtained from the optimized solutions perform better. The last row shows the value of the objective function (here strain energy) for each of the bead patterns.

|             Manual Solution 1              |             Manual Solution 2              |             Manual Solution 3              |             Manual Solution 4              |             Manual Solution 5              |          Optimized Solution           |
| :----------------------------------------: | :----------------------------------------: | :----------------------------------------: | :----------------------------------------: | :----------------------------------------: | :-----------------------------------: |
| <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_1_2D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_2_2D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_3_2D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_4_2D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_5_2D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Opt_2D.png"> |
| <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_1_3D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_2_3D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_3_3D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_4_3D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Manual_5_3D.png"> | <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_Opt_3D.png"> |
|                  2.46e-3                   |                  1.56e-3                   |                  1.17e-3                   |                  1.06e-3                   |                  7.86e-4                   |                4.23e-4                |

### Convergence
The below plots shows the evolution of the objective function over the bead optimization iterations.

<p align="center">
    <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/06_bead_optimization_plate/images/beadOpt_plot.svg" height="650">
</p>


## Source: 
[https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/06_bead_optimization_plate](https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/06_bead_optimization_plate)
