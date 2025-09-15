---
title: Control
keywords: 
tags: [control, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

```Control``` is used to regularize the design variables updates, which are given by the ```Algorithm``` as illustrated in Figure 1. It also updates the properties of the ```model_part```.  This is also used to transform data between physical and control spaces.
<p align="center">
    <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/OptimizationApplication/General/control.png?raw=true" alt="Control data flow"/>
</p>
<p align="center">Figure 1: Control data flow</p>

## Working space

As explained above ```Control``` will work in the control space and physical space.

## Data flow and work flow

```Control``` is responsible for updating the ```model_part```and its properties (i.e. mesh, domain values) with the new design given by the ```Algorithm```. First it will map the control space design (i.e. $$\underline{\hat{\phi}}$$) to physical space (i.e. $$\underline{\phi}$$) using the ```Filtering```. Then these physical space values are used to update the mesh or the domain values.

```Control``` is also responsible for mapping physical space gradients (i.e. $$\frac{dJ_1}{d\underline{\phi}}$$) to control space (i.e. $$\frac{dJ_1}{d\underline{\hat{\phi}}}$$) using the ```Filtering``` method as illustrated in the Figure 1.

## Notes

1. A control should only work in one domain (or ```Kratos::ModelPart```). If more than one domain is required control (or ```Kratos::ModelPart```) then, they should made to one domain using union utilities.

## Source files
* [applications/OptimizationApplication/python_scripts/controls/control.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/controls/control.py)
* [Doxygen](TODO) TODO

