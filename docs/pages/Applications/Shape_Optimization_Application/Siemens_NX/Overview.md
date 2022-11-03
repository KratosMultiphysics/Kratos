---
title: Overview
keywords: 
tags: [NXOverview.md]
sidebar: shape_optimization_application
summary: 
---

The capabilities of the **KratosMultiphysics** **ShapeOptimizationApplication** is exposed in **Siemens NX** through this plugin.

<p align="center">
    <img src="images/plugin_ribbon.png" alt="Plugin"/>
</p>
<p align="center">Figure 1: Plugin ribbon</p>

Figure 1 illustrates the ribbon available in this plugin which is exposing different functionalities of the **ShapeOptimizationApplication**.

## Objectives

This section defines the objectives of the optimization problem. Available objectives are as followings.

1. [Mass](Responses/mass.html)
2. [Strain energy](Responses/strain_energy.html)
3. [Kreisselmeier aggregated stress](Responses/kreisselmeier_aggregated_stress.html)

## Constraints

Each of the above objectives can also be used as constraints in an optimization problem. The specific settings for constraints are explained in [constraints](Objectives_and_constraints.html) section. The list of constraints:

1. [Mass](Responses/mass.html)
2. [Strain energy](Responses/strain_energy.html)
3. [Kreisselmeier aggregated stress](Responses/kreisselmeier_aggregated_stress.html)

## Algorithms

The implemented optimization procedure supports optimization of non-linear objectives under non-linear constraints as well. The options available in thse algorithms are explained in [optimization settings](Optimization_settings.html). Following algorithms may be used:

1. [Steepest descent](Optimization_settings.html#steepest-descent)
2. [Penalized projection](Optimization_settings.html#penalized-projection)

## Design variables

This section determines which are the design variables which are allowed to be controlled by the optimization algorithm to obtain optimum design surfaces. The settings specific to this section can be found in [design variables](Design_variables.html) section.