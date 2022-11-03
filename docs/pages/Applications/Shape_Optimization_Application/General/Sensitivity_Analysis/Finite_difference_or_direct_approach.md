---
title: Finite difference or direct approach
keywords: 
tags: [Finite_difference_or_direct_approach.md]
sidebar: shape_optimization_application
summary: 
---
To calculate the sensitivity of an objective function F to a modification of the design parameters **s**, we start from the unconstrained optimization problem

![f1]

where **u** are the state variables that depend on **s** and **S** is the set of state equations of the system. Applying the chain rule of differentiation to F and S and operating, this yields:

![f2]

where (·)<sub>a</sub> indicates the partial derivative with respect to a. The evaluation of these partial derivatives is problem-specific. In most cases, they are calculated following a finite differences scheme where $$\Delta x_i$$ is the perturbation size:

![f3]

This method is extremely inefficient when the number of degrees of freedom greatly exceeds the number of functions to be evaluated, as is usually the case in shape optimization problems. Therefore, this approach is not recommended except for simple cases.

### References
- Bletzinger, K.-U. (2017). Shape Optimization. In Encyclopedia of Computational Mechanics Second Edition (eds E. Stein, R. Borst and T.J.R. Hughes). https://doi.org/10.1002/9781119176817.ecm2109

[f1]: https://latex.codecogs.com/png.image?\inline&space;\dpi{110}\bg{white}\begin{align*}&space;f=f(\mathbf{s},\mathbf{u}(\mathbf{s})),&space;\\&space;&space;\underline{R}=\underline{R}(\mathbf{s},\mathbf{u}(\mathbf{s}))&space;\end{align*}
[f2]:https://latex.codecogs.com/png.image?\inline&space;\dpi{110}\bg{white}\frac{dF}{d\mathbf{s}}=f_s-f_u\underline{R}_u^{-1}\underline{R}_s
[f3]: https://latex.codecogs.com/png.image?\inline&space;\dpi{110}\bg{white}\dfrac{df}{dx_i}=\dfrac{f(x_i&plus;\Delta&space;x_i)-f(x_i)}{\Delta&space;x_i}