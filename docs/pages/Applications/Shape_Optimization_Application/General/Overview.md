---
title: Overview
keywords: 
tags: [overview.md]
sidebar: shape_optimization_application
summary: 
---

The shape optimization application consist of response functions and algorithms required to obtain optimum design adhering to given constraints and minimizing the given objectives.

It consist of the following technologies:

## [Vertex morphing](../Technologies/Vertex_morphing.html)

This technology is used to minimize the noise in the computed sensitivity fields specially when the objective or the constraint is highly non-linear.

## Optimization algorithms

Following algorithms are used to obtain optimum designs for given objectives satisfying the given constraints.

1. [Steepest descent](../Technologies/Algorithms/steepest_descent.html)
2. [Penalized projection](../Technologies/Algorithms/penalized_projection.html)

## Siemens NX Plugin

The capabilities of $KratosMultiphysics$ $ShapeOptimization$ application is exposed to $Siemens NX$ via a [plugin](../Siemens%20NX/Overview.html)

## Examples

There are two types of example.

1. [Showcase examples](../Examples)
2. [Guided examples](../Siemens_NX/Guided_beam_example/Primal_problem_construction.html)