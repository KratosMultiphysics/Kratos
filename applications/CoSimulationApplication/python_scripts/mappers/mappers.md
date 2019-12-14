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

![\boldsymbol{\Phi} \cdot \boldsymbol{c} = \boldsymbol{\Phi}_{to}](https://render.githubusercontent.com/render/math?math=%5Cboldsymbol%7B%5CPhi%7D%20%5Ccdot%20%5Cboldsymbol%7Bc%7D%20%3D%20%5Cboldsymbol%7B%5CPhi%7D_%7Bto%7D)



[//]: # (MarkDown cheat sheet: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#tables)

[//]: # (generate LaTeX with https://alexanderrodin.com/github-latex-markdown/)