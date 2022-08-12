---
title: Multi objective optimization
keywords: 
tags: [Multi_objective_optimization.md]
sidebar: shape_optimization_application
summary: 
---
Multi-objective or multi-criteria optimization takes place when more than one objective function is characterized in the optimization problem. In general, these objectives are conflicting, requiring to find a compromise between the possible solutions.

## Compromise function approach

A first approach would be combining the different objective functions in a weighted sum:

![f1]

The choice of the weights w<sub>i</sub> determine the importance of each objective function f<sub>i</sub> and determines the optimal solution. The constant C is introduced to shift the compromise function, but does not affect the optimal solution.

## Hierarchical approach
A second option is choosing one of the objective functions as the main one and rewrite the rest as constrains. Limits can be introduced in those constraints to ensure that the obtained values lie on a certain range.


## Pareto optimization
Pareto's optimization strategy follows a systematic approach to find the compromise solution. A solution is considered Pareto optimal if, for any modification of the design variables, at least one of the objectives functions becomes worse. The so-called Pareto front groups all the points that are considered Pareto optimal for a given multi-objective problem. 

The construction of the Pareto front is usually challenging and computationally expensive, as it requires the evaluation of many problem configurations. Then, those cases that are not dominated by any other are chosen as the Pareto front. It is said that a case is dominated if there exist another configuration that provides a gain in some objective functions and none of the others lose. Finally, a compromise point is chosen among those in the Pareto front.

To choose the compromise point, the ideal optimum made up of the individual optimums of the objective function is taken as reference. The quality of a solution f<sub>i</sub>(x) can be defined by its normalized distance d<sub>i</sub> from the individual optimum f<sub>i</sub><sup>*</sup>(x):
 
![f2]

The best compromise can be defined by the minimum of the largest deviation (min-max optimum) or by the minimum distance from the ideal optimum (Euclidean norm optimum).


# Kratos implementation

Multi-objective optimization is not directly implemented in Kratos. Nevertheless, the simulations results obtained in Kratos can be used in combination with the previous techniques. 



[f1]: https://latex.codecogs.com/png.image?\inline&space;\dpi{110}\bg{white}F(\mathbf{x})&space;=&space;\sum_{i=1}^{p}w_if_i(\mathbf{x})&plus;C
[f2]: https://latex.codecogs.com/png.image?\dpi{110}\bg{white}d_i(\mathbf{x})=\dfrac{f_i(\mathbf{x})-f_i^*}{f_i^*}