---
title: Algorithm
keywords: 
tags: [algorithm, optimization]
sidebar: optimization_application
summary: 
---

## Introduction

Figure 1 illustrates how an algorithm operates when ```OptimizationAnalysis``` executes ```Algorithm::Solve``` method. Overview of the existing algorithms can be found [here](../Algorithms/Overview.html).

<p align="center">
    <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/OptimizationApplication/General/algorithm.png?raw=true" alt="Algorithm"/>
</p>
<p align="center">Figure 1: Algorithm</p>

## Working space.

As illustrated in the Figure 1, ```Algorithm``` purely works in the control space. Hence, all the inputs to the ```Algorithm``` should be in the control space. Therefore, all the outputs from it also will be in the control space.

## Supporting components

When an ```Algorithm``` is constructed by the ```OptimizationAnalysis```,
1. It will create the ```MasterControl``` which links all the ```Control```s used by the ```Algorithm```.
2. It will create all the necessary ```ResponseRoutine```s for the algorithm. It will not create any ```ResponseFunction```s, it will merely link a ```ResponseRoutine``` with a ```ResponseFunction```.
3. Thereafter, it will create all the other necessary components.

## Data flow and work flow

```Algorithm``` is developed basically to perform minimization. Therefore, it is the duty of the ```ResponseRoutine``` to standardize the objective and constraint values (i.e. $$\tilde{J}_1$$ and $$\tilde{J}_2$$ and the gradients (i.e. $$\frac{d\tilde{J}_1}{d\underline{\hat{\phi}}}$$, $$\frac{d\tilde{J}_2}{d\underline{\hat{\phi}}}$$) accordingly and send it to the algorithm which will be based on physical space response values (i.e. $$J_1$$, $$J_2$$) and their gradients (i.e. $$\frac{dJ_1}{d\underline{\phi}}$$, $$\frac{dJ_2}{d\underline{\phi}}$$). Then using the algorithm, it will compute the new design in control space (i.e. $$\hat{\underline{\phi}}$$) which needs to mapped to a design in the physical space (i.e. $$\underline{\phi}$$) by the ```ResponseRoutine```. ```ResponseRoutine``` then should update its design surfaces using ```MasterControl``` for the next iteration computations.


## Notes

1. Since each algorithm is specific, they can have their own specific types of ```ResponseRoutine```s developed. Therefore, the methods being called by an ```Algorithm``` to ```ResponseRoutine``` may differ depending on the type of the ```Algorithm```.
2. ```OptimizationApplication``` supports linking with third party libraries for various algorithms. In that case, one needs to develop the ```Algorithm``` wrapper, ```ResponseRoutine``` and ```MasterControl``` for that third party library.

## Source files
* [applications/OptimizationApplication/python_scripts/algorithms/algorithm.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/algorithms/algorithm.py)
* [Doxygen](TODO) TODO


