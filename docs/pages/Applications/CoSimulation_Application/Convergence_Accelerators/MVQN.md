---
title: MVQN
keywords: 
tags: [MVQN.md]
sidebar: cosimulation_application
summary: 
---
## Overview
The Multi Vector Quasi Newton is a formulation based on Newton-Raphson method. This formulation approximate the Jacobian, based on the residual and ∆𝑢 of the previous iterations.


The Taylor expansion of the residual is expressed as:
<p align="center">$$ r^{i+1}\approx r^i+\frac{\partial r^i}{\partial u}\Delta u^{i+1}$$</p> 

where $\frac{\partial r^i}{\partial u}$ is the jacobian $J^i$.
Based on the equation above, the residual will be zero when:
<p align="center">$$ \frac{\partial r^i}{\partial u}\Delta u^{i+1}=-r^i$$</p>

The jacobian is aproximated with the following equation:
<p align="center">$$ J^{i+1,t+1}\approx J^t+(W^i-J^tV^i)((V^i)^TV^i)^{-1}(V^i)^T$$</p>
where V is a collection of residual vectors from the previous iterations and W is a collection of $\Delta u$ vectors from the previous iterations.

solving for $\Delta u^{i+1}$, we can calculate the next iteration:
<p align="center">$$ u^{i+1}=u^i+\Delta u^{i+1}$$</p>


## Parameters:
- **horizon**: Maximum number of vectors to be stored in each time step.
- **alpha**: Relaxation factor for computing the update, when no vectors available.

## References
- [A.E.J. Bogaers et al. "Quasi-Newton methods for implicit black-box FSI coupling", Computational methods in applied mechanics and engineering. 279(2014) 113-132.](https://www.researchgate.net/publication/263858221_Quasi-Newton_methods_for_implicit_black-box_FSI_coupling)