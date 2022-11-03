---
title: Strain energy
keywords:
tags: [Strain_energy.md]
sidebar: shape_optimization_application
summary:
---

## Introduction

This computes the summation of strain energy from each element as the response value and its gradient.

## Formulation

Following formulation is used to compute the summation of strain energy from each element in the chosen model part where $$\underline{u}$$ is the displacement vector for the element and $$\mathbf{K}$$ is the stiffness matrix of the element.
<p align="center">$$ J   = \frac{1}{2}\sum_{\Omega} \underline{u}^T \mathbf{K} \underline{u}  $$</p>

Location: ["applications/StructuralMechanicsApplication/python_scripts/structural_response.py"](https://github.com/KratosMultiphysics/Kratos/blob/shapeopt/kreisselmeier_aggregation/applications/StructuralMechanicsApplication/python_scripts/structural_response.py)
