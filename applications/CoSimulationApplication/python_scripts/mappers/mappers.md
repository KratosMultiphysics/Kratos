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
`direction`|str|coordinate direction, options: `'X'`, `'Y'`, `'Z'`
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build


#### MapperLinear2D

JSON setting|type|description
------:|:----:|-----------
`direction_1`|str|first coordinate direction, options: `'X'`, `'Y'`, `'Z'`
`direction_2`|str|second coordinate direction, options: `'X'`, `'Y'`, `'Z'`
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build

#### MapperLinear3D

JSON setting|type|description
------:|:----:|-----------
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build
`parallel`|bool|if `true`, use `multiprocessing` module to distribute calculation of coefficients over all cores of current node


#### MapperRadialBasis2D

With φ(r) a radial basis function defined as  

> φ(r) = (1 − r)<sup>4</sup> (1 + 4r) for 0 ≤ r < 1  
> φ(r) = 0 for 1 ≤ r

a function _f_(**x**) can be approximated as the weighted sum of _n_ shifted radial basis functions:

> _f_(**x**) = Σ<sub>j</sub> α<sub>j</sub> φ(||**x** − **x**<sub>j</sub>||)

To determine the coefficients, we require that the exact function value is returned at the _n_ test points.
This gives us _n_ equations

> _f_(**x**<sub>i</sub>) = **f**<sub>i</sub> 
= Σ<sub>j</sub> α<sub>j</sub> φ(||**x**<sub>i</sub> − **x**<sub>j</sub>||)

which can be written in matrix form as

> **f** = **Φ** · **α**

with **f**, **α** ∈ R<sup>n×1</sup> and **Φ** ∈ R<sup>n×n</sup>. 
The weights-vector **α** can be extracted by solving this system:

> **α** =  **Φ**<sup>-1</sup> · **f**





![\boldsymbol{\Phi} \cdot \boldsymbol{c} = \boldsymbol{\Phi}_{to}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D%20%5Ccdot%20%5Cboldsymbol%7Bc%7D%20%3D%20%5Cboldsymbol%7B%5CPhi%7D_%7Bto%7D)



[//]: # (MarkDown cheat sheet: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#tables)

[//]: # (generate LaTeX eqns: https://alexanderrodin.com/github-latex-markdown/)

[//]: # (HTML math symbols: http://www.unics.uni-hannover.de/nhtcapri/mathematics.html)
[//]: # (more: http://www.alanflavell.org.uk/unicode/unidata22.html)

[//]: # (Greek lower: αβγδεζηϑθικλμνξοπρστυφϕχψω)
[//]: # (Greek upper: ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ)
[//]: # (super, sub: <sup></sup>, <sub></sub> )
[//]: # (operators: + - − · / × √ ∘ ∗)
[//]: # (other: ∂ Δ	∑ ≤ ≥ ∈ )