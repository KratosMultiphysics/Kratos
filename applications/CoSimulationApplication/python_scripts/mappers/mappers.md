# Mappers



## General concepts


#### Interpolators and Transformers

#### Calling the mapper


## Overview of available mappers

#### MapperInterface

Special mapper-class: takes two `CoSimulationInterface` objects, 
and maps the `ModelPart` objects to each other in order of appearance. 


#### MapperLinear1D

JSON setting|type|description
------:|:----:|-----------
`direction`|str|coordinate direction, options: `"X"`, `"Y"`, `"Z"`
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build


#### MapperLinear2D

JSON setting|type|description
------:|:----:|-----------
`direction_1`|str|first coordinate direction, options: `"X"`, `"Y"`, `"Z"`
`direction_2`|str|second coordinate direction, options: `"X"`, `"Y"`, `"Z"`
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build

#### MapperLinear3D

JSON setting|type|description
----------:|:---:|-----------
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build
`parallel`|bool|if `true`, use `multiprocessing` module to distribute calculation of coefficients over all cores of current node


#### MapperRadialBasis

JSON setting|type|description
----------:|:---:|-----------
`directions`|str/list|one or more directions to interpolate, given as a list consisting of `"X"`, `"Y"`, `"Z"`; for 1D, a string is also accepted
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build



With φ(r) a radial basis function defined as  

> φ(r) = (1 − r)<sup>4</sup> (1 + 4r) for 0 ≤ r < 1  
> φ(r) = 0 for 1 ≤ r

![\phi(r)=\begin{cases}(1-r)^4(1+4r)&r\lt1\\0&r\geq1\end{cases}](https://render.githubusercontent.com/render/math?math=%5Cphi(r)%3D%5Cbegin%7Bcases%7D(1-r)%5E4(1%2B4r)%26r%5Clt1%5C%5C0%26r%5Cgeq1%5Cend%7Bcases%7D)

a function _f_(**x**) can be approximated as the weighted sum of _n_ shifted radial basis functions:

> _f_(**x**) = Σ<sub>j</sub> α<sub>j</sub> φ(||**x** − **x**<sub>j</sub>||)

![f(\boldsymbol{x})=\sum_j^n\alpha_j\,\phi(||\boldsymbol{x}-\boldsymbol{x}_j||) ](https://render.githubusercontent.com/render/math?math=f(%5Cboldsymbol%7Bx%7D)%3D%5Csum_j%5En%5Calpha_j%5C%2C%5Cphi(%7C%7C%5Cboldsymbol%7Bx%7D-%5Cboldsymbol%7Bx%7D_j%7C%7C)%20)


To determine the coefficients, we require that the exact function value is returned at the _n_ test points.
This gives us _n_ equations

> _f_(**x**<sub>i</sub>) = **f**<sub>i</sub> 
= Σ<sub>j</sub> α<sub>j</sub> φ(||**x**<sub>i</sub> − **x**<sub>j</sub>||)

![f(\boldsymbol{x}_i)=\sum_j^n\alpha_j\,\phi(||\boldsymbol{x}_i-\boldsymbol{x}_j||) ](https://render.githubusercontent.com/render/math?math=f(%5Cboldsymbol%7Bx%7D_i)%3D%5Csum_j%5En%5Calpha_j%5C%2C%5Cphi(%7C%7C%5Cboldsymbol%7Bx%7D_i-%5Cboldsymbol%7Bx%7D_j%7C%7C)%20)


which can be written in matrix form as

> **f** = **Φ** · **α**

![\boldsymbol{f}=\boldsymbol{\Phi}\cdot\boldsymbol{\alpha}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7Bf%7D%3D%5Cboldsymbol%7B%5CPhi%7D%5Ccdot%5Cboldsymbol%7B%5Calpha%7D)


with **f**, **α** ∈ R<sup>n×1</sup> and **Φ** ∈ R<sup>n×n</sup>. 
The weights-vector **α** can be extracted by solving this system:

> **α** =  **Φ**<sup>-1</sup> · **f**

![\boldsymbol{\alpha}=\boldsymbol{\Phi}^{-1}\cdot\boldsymbol{f}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5Calpha%7D%3D%5Cboldsymbol%7B%5CPhi%7D%5E%7B-1%7D%5Ccdot%5Cboldsymbol%7Bf%7D)

However, the function value vector **f** is not known in advance: it contains the values of the `Variable` that will be interpolated. 
What we do want, is the weights (coefficients) with which every value in the _n_ nearest from-points must be multiplied, to get the value in the to-point. 

...


![\boldsymbol{\Phi} \cdot \boldsymbol{c} = \boldsymbol{\Phi}_{to}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D%20%5Ccdot%20%5Cboldsymbol%7Bc%7D%20%3D%20%5Cboldsymbol%7B%5CPhi%7D_%7Bto%7D)



[//]: # (MarkDown cheat sheet: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#tables)

[//]: # (render LaTeX eqn as image: https://alexanderrodin.com/github-latex-markdown/)

[//]: # (HTML math symbols: http://www.unics.uni-hannover.de/nhtcapri/mathematics.html)
[//]: # (more: http://www.alanflavell.org.uk/unicode/unidata22.html)

[//]: # (Greek lower: αβγδεζηϑθικλμνξοπρστυφϕχψω)
[//]: # (Greek upper: ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ)
[//]: # (super, sub: <sup></sup>, <sub></sub> )
[//]: # (operators: + - − · / × √ ∘ ∗)
[//]: # (other: ∂ Δ	∑ ≤ ≥ ∈ )