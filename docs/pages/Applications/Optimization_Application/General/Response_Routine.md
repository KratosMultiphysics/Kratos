---
title: Response Routine
keywords: 
tags: [response routine, response function, optimization]
sidebar: optimization_application
summary: 
---
## Introduction

Figure 1 illustrates how a ```ResponseRoutine``` operates.

<p align="center">
    <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/OptimizationApplication/General/response_routine.png?raw=true" alt="Response routine"/>
</p>
<p align="center">Figure 1: Response routine</p>

## Working space.

As illustrated in the Figure 1, ```ResponseRoutine``` will work in the control space and physical space. It uses ```MasterControl``` to transform control space to physical space and wise versa.

## Supporting components

Each ```ResponseRoutine``` will not create its own ```ResponseFunction```s. It should use ```ResponseFunction```s given by the ```OptimizationProblem``` data container. Once assigned, ```ResponseRoutine``` has the authority on how to control the ```ResponseFunction``` (such as calling methods for ```ResponseFunction```).

**```ResponseRoutine``` does not own the ```MasterControl```. It should merely use the shared ```MasterControl```.**

## Data flow and work flow

```ResponseRoutine``` will receive new design in the control space (i.e. $$\underline{\hat{\phi}}$$). Then it should use the ```MasterControl``` to convert the new design to physical space (i.e. $$\underline{\phi}$$) and update the mesh. Thereafter, it will request the new objective/constraint values (i.e. $$J_1$$) and their gradients in the physical space (i.e. $$\frac{dJ_1}{d\underline{\phi}}$$). Thereafter, it will convert the physical space gradients to control space gradients (i.e. $$\frac{dJ_1}{d\underline{\hat{\phi}}}$$). Finally it will return the standardized objective/constraint values (i.e. $$\tilde{J}_1$$) and their gradients (i.e. $$\frac{d\tilde{J}_1}{d\underline{\hat{\phi}}}$$)

## Notes
* One of the major duties of the ```ResponseRoutine``` is to request the physical control variables
used in the ```MasterControl``` for each domain, and find out what gradients should be requested from the ```ResponseFunction```s.

## Source files
* [applications/OptimizationApplication/python_scripts/responses/response_routine.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/responses/response_routine.py)
* [Doxygen](TODO) TODO