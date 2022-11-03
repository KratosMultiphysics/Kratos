---
title: Damping
keywords: 
tags: [Damping.md]
sidebar: shape_optimization_application
summary: 
---

## Introduction

Damping of sensitivity gradients are used to reduce or restrict the shape update getting changed due to objective gradient anr/or active constraint gradients.

### Damping functions

Following damping functions are supported:

#### Gaussian filter function ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/custom_utilities/filter_function.cpp#L34))

<p align="center">$$ A\left(\mathbf{x},\mathbf{x_0},r\right)  = \max\left\lbrace 0.0, e^{\frac{-9\left|\mathbf{x}-\mathbf{x_0}\right|^2}{2r^2}}\right\rbrace$$</p>

#### Linear filter function ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/custom_utilities/filter_function.cpp#L38))
<p align="center">$$ A\left(x,x_0,r\right)  = \max\left\lbrace 0.0, \frac{r-\left|\mathbf{x}-\mathbf{x_0}\right|}{r}\right\rbrace$$</p>

#### Constant filter function ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/custom_utilities/filter_function.cpp#L42))
<p align="center">$$ A\left(x,x_0,r\right)   = 1.0$$</p>

#### Cosine filter function ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/custom_utilities/filter_function.cpp#L46))
<p align="center">$$ A\left(x,x_0,r\right)   = \max\left\lbrace 0.0, 1-0.5\left(1-\cos\left(\pi\frac{\left|\mathbf{x}-\mathbf{x_0}\right|}{r}\right)\right)\right\rbrace$$</p>
<!-- [](double radius, double distance) {return std::max(0.0, );}; -->

#### Quartic filter function ([source](https://github.com/KratosMultiphysics/Kratos/blob/0048ec0790af5b356039ee4829d78ff0deb2d640/applications/ShapeOptimizationApplication/custom_utilities/filter_function.cpp#L50))
<p align="center">$$ A\left(x,x_0,r\right)   = \max\left\lbrace 0.0, \left(\frac{\left|\mathbf{x}-\mathbf{x_0}\right|-r}{r}\right)^4\right\rbrace$$</p>

## Damping methodology

Damping factors (i.e. $$\beta_{i, damp}$$) are computed for each $$i^{th}$$ node in the damping regions as follows. First neighbours of each node in the damping region is found using a KDTree search.
<p align="center">$$ \beta_{i, damp} = \min_{x_j \in \Omega_{i, neighbours}} 1.0 - A(x_i, x_j, r) $$</p>

Thereafter, the sensitivity of the node is multiplied by the damping factor to compute the damped sensitivities as shown in the following equation:
<p align="center">$$ \left(\frac{dh}{d\underline{s}}\right)_{i, damped} = \beta_{i, damp}\left(\frac{dh}{d\underline{s}}\right)_i $$</p>


## Source

Location: ["applications/ShapeOptimizationApplication/custom_utilities/damping/damping_utilities.h"](https://github.com/KratosMultiphysics/Kratos/blob/shapeopt/kreisselmeier_aggregation/applications/ShapeOptimizationApplication/custom_utilities/damping/damping_utilities.h)
