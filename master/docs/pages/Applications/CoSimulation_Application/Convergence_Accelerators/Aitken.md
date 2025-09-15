---
title: Aitken
keywords: 
tags: [Aitken.md]
sidebar: cosimulation_application
summary: 
---
Coupling convergence can be improved by applying a relaxation parameter. Ideally, it has to be small enough to keep the iteration from diverging, but as large as possible to minimize the amount of iteration. 
<p align="center">$$ u_1^{i + 1} = \hat u_2^{i + 1} = \alpha u_2^{i + 1} + (1 - \alpha )u_2^i$$</p>
To maximize the efficiency, Aitken's formulation dynamically calculates a suitable relaxation parameter for each iteration by taking into account the result from the previous iterations.
<p align="center">$$ {\alpha ^{i + 1}} = {\alpha ^i}\frac{(r^i)^T(r^{i+1}-r^i)}{||r^{i+1}-r^i||}$$</p>
<p align="center">$$ r^{i+1}=u_2^{i+1}-u_1^{i+1}$$</p>

The default parameters specific to this solver are:
```json
    "init_alpha"     :  0.1,
    "init_alpha_max" :  0.45,
    "alpha_max"      :  2.0,
    "alpha_min"      : -2.0
```
- **init_alpha**: Relaxation factor in the first time step.
- **init_alpha_max**: Maximum relaxation factor for the first iteration in each time step.
- **alpha_max**: Upper bound for the dynamic relaxation factor. 
- **alpha_min**: Lower bound for the dynamic relaxation factor.

Reference: [Ulrich Küttler et al., "Fixed-point fluid–structure interaction solvers with dynamic relaxation"](https://link.springer.com/article/10.1007/s00466-008-0255-5#article-info)