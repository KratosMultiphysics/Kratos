---
title: Execution Policy
keywords: 
tags: [execution policy, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

Execution policies are used to execute an analysis which will be used to compute response values
and/or gradients. Overview of available ```ExecutionPolicy```s can be found [here](../Execution_Policies/Overview.html).

These analysis may be:
* Kratos static analysis
* Kratos dynamic analysis
* External analysis

Depending on the analysis, you may have to develop your won ```ExecutionPolicy``` to be used with optimization analysis.

All the information about new design updates are passed via ```Kratos::ModelPart```, hence no additional information will be passed to this except for ```OptimizationProblem``` data container.

## Source files
* [applications/OptimizationApplication/python_scripts/execution_policies/execution_policy.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/execution_policies/execution_policy.py)
* [Doxygen](TODO) TODO
