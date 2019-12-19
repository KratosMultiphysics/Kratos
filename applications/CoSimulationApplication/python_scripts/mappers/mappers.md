# Mappers


## General concepts


#### Hierarchy of mapping-related objects

CoCoNuT interacts with the mappers through the `SolverWrapperMapped` object: this wrapper behaves like every other `SolverWrapper` as far as the other components are concerned. 
It contains 3 main components: a `Mapper` for the input, a real `SolverWrapper` and a `Mapper` for the output. The mappers are initialized through the `SetInterfaceInput` and `SetInterfaceOutput` methods respectively, by providing them with the `CoSimulationInterface` objects that will be respectively the input and output of the `SolverWrapperMapped` object.

The two mappers in the `SolverWrapperMapped` object are also of a special type: they work on the level of `CoSimulationInterface` objects. They are some sort of mapper-wrapper around the actual mappers which work on `ModelPart` level.
Currently only one such mapper is available, aptly called `MapperInterface`.

At the lowest level, mappers interpolate historical variables between two `ModelPart` objects, based on the coordinates of the nodes.
They can be chained together in a `MapperCombined` object, creating in fact another layer of mapping. So many layers! Like an onion!


#### Interpolators and transformers

The `ModelPart`-level mappers have two important methods.

First, an `Initialize` method in which one-time expensive operations are performed, mostly nearest-neighbour searches and calculation of interpolation coefficients. This initialization is done based on the original coordinates `X0`, `Y0` and `Z0`. 

Second, a `__call__` method which is used for the actual mapping. It takes two tuples as arguments, containing the `ModelPart` and `Variable` which are used in the interpolation. This method returns nothing: the interpolation is done in-place in the `ModelPart` objects.

There are two types of `ModelPart`-level mappers: interpolators and transformers. They can be distinguished by their boolean `interpolator` attribute (see `__init__`). 

For interpolators, `Initialize` gets two `ModelPart` objects (dubbed _from_ and _to_), and returns nothing. These mappers do real interpolation, examples are `MapperNearest` and `MapperRadialBasis`.

For transformers, `Initialize` gets only one `ModelPart` (from either the _from_ or _to_ side, depending on the transformation), and returns the other `ModelPart`. An example is the `MapperPermutation` transformer, which exchanges the coordinates as specified. 

A transformer can never be used by itself, it must always be combined with an interpolator: the reason is that interpolators use information coming from two sides, which is exactly what the `SolverWrapperMapped` and `MapperInterface` objects want. To chain together multiple mappers, the `MapperCombined` is used: it contains always 1 interpolator and 0 or more transformers, on either side of the interpolator.


## Overview of available mappers

> TODO: give some details/explanation about every implemented mapper. Shorten comments in the actual Python code.

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


TO DO: write this for only 1 point, because nearest neighbours differ!

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


with **f**, **α** ∈ R<sup>n×1</sup>, 
![\boldsymbol{f},\boldsymbol{\alpha}\in\mathbb{R}^{n\times1}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7Bf%7D%2C%5Cboldsymbol%7B%5Calpha%7D%5Cin%5Cmathbb%7BR%7D%5E%7Bn%5Ctimes1%7D)
and **Φ** ∈ R<sup>n×n</sup>
![\boldsymbol{\Phi}\in\mathbb{R}^{n\times n}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D%5Cin%5Cmathbb%7BR%7D%5E%7Bn%5Ctimes%20n%7D)
. The weights-vector **α** can be extracted by solving this system:

> **α** =  **Φ**<sup>-1</sup> · **f**

![\boldsymbol{\alpha}=\boldsymbol{\Phi}^{-1}\cdot\boldsymbol{f}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5Calpha%7D%3D%5Cboldsymbol%7B%5CPhi%7D%5E%7B-1%7D%5Ccdot%5Cboldsymbol%7Bf%7D)

However, the function value vector **f** is not known in advance: it contains the values of the `Variable` that will be interpolated. 
What we do want, is the weights (coefficients) with which every value in the _n_ 
![n](https://render.githubusercontent.com/render/math?math=n) 
nearest from-points must be multiplied, to get the value in the to-point. 
The target value _f_<sub>to</sub> would be calculated as

![f(\boldsymbol{x}_{to})=\sum_j^n\alpha_j\,\phi(||\boldsymbol{x}_{to}-\boldsymbol{x}_j||)=\boldsymbol{\Phi}^T_{to}\cdot\boldsymbol{\alpha}=\boldsymbol{\Phi}^T_{to}\cdot\boldsymbol{\Phi}^{-1}\cdot\boldsymbol{f}=\boldsymbol{c}^T\cdot\boldsymbol{f}](https://render.githubusercontent.com/render/math?math=f(%5Cboldsymbol%7Bx%7D_%7Bto%7D)%3D%5Csum_j%5En%5Calpha_j%5C%2C%5Cphi(%7C%7C%5Cboldsymbol%7Bx%7D_%7Bto%7D-%5Cboldsymbol%7Bx%7D_j%7C%7C)%3D%5Cboldsymbol%7B%5CPhi%7D%5ET_%7Bto%7D%5Ccdot%5Cboldsymbol%7B%5Calpha%7D%3D%5Cboldsymbol%7B%5CPhi%7D%5ET_%7Bto%7D%5Ccdot%5Cboldsymbol%7B%5CPhi%7D%5E%7B-1%7D%5Ccdot%5Cboldsymbol%7Bf%7D%3D%5Cboldsymbol%7Bc%7D%5ET%5Ccdot%5Cboldsymbol%7Bf%7D)

where **c** is calculated by solving the system

![\boldsymbol{\Phi}\cdot\boldsymbol{c}=\boldsymbol{\Phi}_{to}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D%5Ccdot%5Cboldsymbol%7Bc%7D%3D%5Cboldsymbol%7B%5CPhi%7D_%7Bto%7D).

As every to-point has different nearest neighbours in the from-points, the coefficient vector **c** must be calculated for each to-point independently. The matrices 
![\boldsymbol{\Phi}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D) 
and 
![\boldsymbol{\Phi}_{to}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D_%7Bto%7D) 
must also be recalculated for every to-point.




[//]: # (MarkDown cheat sheet: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#tables)

[//]: # (render LaTeX eqn as image: https://alexanderrodin.com/github-latex-markdown/)

[//]: # (HTML math symbols: http://www.unics.uni-hannover.de/nhtcapri/mathematics.html)
[//]: # (more: http://www.alanflavell.org.uk/unicode/unidata22.html)

[//]: # (Greek lower: αβγδεζηϑθικλμνξοπρστυφϕχψω)
[//]: # (Greek upper: ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ)
[//]: # (super, sub: <sup></sup>, <sub></sub> )
[//]: # (operators: + - − · / × √ ∘ ∗)
[//]: # (other: ∂ Δ	∑ ≤ ≥ ∈ )