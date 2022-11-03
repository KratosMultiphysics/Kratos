---
title: Optimization Settings
keywords: 
tags: [Optimization_settings.md]
sidebar: shape_optimization_application
summary: 
---

## Line search

Line search is used to determine the next shape update along the surface normal sensitivities direction which are computed from the adjoint solution.

Following settings can be found in the `Line search settings` section under the `Optimization settings`.

`Type` contains the methodologies supported in determining next search direction. Currently `manual` is the only methodology supported where surface normal sensitivities are normalized and afterwards they are applied with the manually specified step size (give by `Step size`) to obtain next shape update.

## Algorithms

Mainly there are two algorithms supported in the Siemens NX plugin. They are `Steepest descent` and `Penalized projection`. Following sections explain the settings found in each algorithm.

### Steepest descent

<p align="center">
    <img src="images/algorithm_steepest_descent.png" alt="Steepest descent settings"/>
</p>
<p align="center">Figure 1: Steepest descent settings</p>

Figure 1 illustrates basic settings available in the steepest descent algorithm. Further details of the settings as well as the implementation can be found in the [steepest descent](../Technologies/Algorithms/steepest_descent.html) section.

|Option|Description|
|------|-----------|
|Max. Steps| Maximum number of steps allowed in the optimization procedure|
|Relative Tolerance| Relative tolerance change between two design iterations to check whether optimization procedure has converged|

### Penalized projection

<p align="center">
    <img src="images/algorithm_penalized_projection.png" alt="Penalized projection settings"/>
</p>
<p align="center">Figure 2: Penalized projection settings</p>

Figure 2 illustrates basic settings available in the penalized projection algorithm. Further details of the settings as well as the implementation can be found in the [penalized projection](../Technologies/Algorithms/penalized_projection.html) section.

|Option|Description|
|------|-----------|
|Correction Scaling| $$\alpha_{corscal}^0$$ value (refer [Correction scaling](../Technologies/Algorithms/penalized_projection.html#computing-correction-scaling)) |
|Adaptive correction scaling| If false, then $$\alpha_{corscal}^0$$ is used for all the design iterations. If true then $$\alpha_{corscal}^n$$ is computed (refer [Correction scaling](../Technologies/Algorithms/penalized_projection.html#computing-correction-scaling)) |
|Max. Steps| Maximum number of steps allowed in the optimization procedure|
|Relative Tolerance| Relative tolerance change between two design iterations to check whether optimization procedure has converged|