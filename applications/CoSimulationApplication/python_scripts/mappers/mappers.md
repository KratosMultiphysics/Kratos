# Mappers


## General concepts


#### Hierarchy of mapping-related objects

CoCoNuT interacts with the mappers through the `SolverWrapperMapped` object: this wrapper behaves like every other `SolverWrapper` as far as the other components are concerned. 
It contains 3 main components: a `Mapper` for the input, a real `SolverWrapper` and a `Mapper` for the output. The mappers are initialized through the `SetInterfaceInput` and `SetInterfaceOutput` methods respectively, by providing them with the `CoSimulationInterface` objects that will be respectively the input and output of the `SolverWrapperMapped` object.

The two mappers in the `SolverWrapperMapped` object are also of a special type: they work on the level of `CoSimulationInterface` objects. They are some sort of mapper-wrapper around the actual mappers which work on `ModelPart` level.
Currently only one such mapper is available, aptly called `MapperInterface`.

At the lowest level, mappers interpolate historical variables between two `ModelPart` objects, based on the coordinates of the nodes. Interpolation is always done from the _from_-`ModelPart` to the _to_-`ModelPart`.
These mappers can be chained together in a `MapperCombined` object, creating in fact another layer of mapping. So many layers! Like an onion!


#### Interpolators and transformers

The `ModelPart`-level mappers have two main methods: `Initialize` and `__call__`. 

The `Initialize` method performs one-time expensive operations, namely nearest-neighbour search and calculation of the interpolation coefficients. The initialization is done based on the original coordinates `X0`, `Y0` and `Z0`.

The `__call__` method is used for the actual mapping. It takes two tuples as arguments (_from_ and _to_ respectively), each tuple containing the `ModelPart` and the `Variable` to be used in the interpolation. This method returns nothing: the interpolation is done in-place in the `ModelPart` objects.

There are two types of `ModelPart`-level mappers: interpolators and transformers. They can be distinguished by their boolean `interpolator` attribute (see `__init__`). 

For interpolators, `Initialize` gets two `ModelPart` objects (_from_ and _to_), and returns nothing. These mappers do the real interpolation. Currently `MapperNearest`, `MapperLinear` and `MapperRadialBasis` are available.

For transformers, `Initialize` gets only one `ModelPart` (from either the _from_ or _to_ side, depending on the transformation), and returns the other `ModelPart`. An example is the `MapperPermutation` transformer, which exchanges the coordinates as specified. 

A transformer can never be used by itself, it must always be combined with an interpolator: the reason is that interpolators use information coming from two sides, which is exactly what the `SolverWrapperMapped` and `MapperInterface` objects want. To chain together multiple mappers, the `MapperCombined` is used: it contains always 1 interpolator and 0 or more transformers, on either side of the interpolator.


## Overview of available mappers


#### MapperInterface

Special mapper-class that maps on the level of `CoSimulationInterface` objects. 
It takes two `CoSimulationInterface` objects, 
and maps the `ModelPart` objects to each other in order of appearance, all using the same `ModelPart` mapper.

To use different interpolation for the different `ModelPart` objects or even for different historical variables, a new `CoSimulationInterface` mapper must be written. 

JSON setting|type|description
------:|:----:|-----------
`type`|str|`ModelPart` mapper to be used
`settings`|dict|all the settings for the `ModelPart` mapper specified in `type`



#### MapperInterpolator

Base-class for all interpolators (currently `MapperNearest`, `MapperLinear` and `MapperRadialBasis`). 

JSON setting|type|description
------:|:----:|-----------
`directions`|list|list of coordinate directions, maximum three entries, may contain `"X"`, `"Y"`, `"Z"`
`balanced_tree`|bool|if `true`, create balanced `cKDTree`, which is more stable, but takes longer to build; set to `true` if the tree is giving problems (which I don't expect)

The `Initialize`-method should be called in all child-classes. It does the following:
- read and store the coordinates from the _from_ and _to_ `ModelPart` objects
- check if the bounding boxes of the _from_ and _to_ `ModelPart` objects are more or less overlapping
- do an efficient nearest neighbour search using `scipy.spatial.cKDTree`
- check if the _from_ `ModelPart` does not contain duplicate nodes (i.e. same coordinates)

The `__call__`-method should not be overridden in the child-classes. It maps historical variables based on neighbours and coefficients determined in `Initialize`. Historical variables of type `Double` and type `Array` can be mapped (the latter is just the application of the former for each vector component).


#### MapperNearest

Child-class of `MapperInterpolator`, does not require additional settings. Does simple nearest-neighbour mapping.


#### MapperLinear

Child-class of `MapperInterpolator`, additional settings:

JSON setting|type|description
------:|:----:|-----------
`parallel`|bool|if `true`, use `multiprocessing` to parallellize loop that calculates coefficients

The kind of linear mapping depends on the number of coordinate directions, as given in the `directions` setting.

**1D** - If the _to_-point lies between the 2 nearest _from_-points, linear interpolation is done. Else, nearest neighbour interpolation is done.

**2D** - The _to_-point is first projected on the line through the 2 nearest _from_-points. If the projected point lies between the _from_-points, linear interpolation is done. Else, nearest neighbour interpolation is done.

**3D** - The _to_-point is first projected on the plane through the 3 nearest _from_-points. If the triangle consinsting of those 3 points is _deprecated_ (colinear points), the 2D-methodology is followed. Else, if the projected point lies inside the triangle, barycentric interpolation is done. If it lies outside the triangle, the 2D-methodology is followed.


#### MapperRadialBasis

Child-class of `MapperInterpolator`, additional settings:

JSON setting|type|description
------:|:----:|-----------
`parallel`|bool|if `true`, use `multiprocessing` to parallellize loop that calculates coefficients

Radial basis function interpolation is relatively straightforward: implementation for 1D, 2D and 3D is exactly the same and can be written in a condensed way using `scipy.spatial.distance`. 

Normal radial basis interpolation is done as follows.
_φ_(_r_) is a radial basis function defined as  

 _φ_(_r_) = (1 − _r_)<sup>4</sup> (1 + 4 _r_) for 0 ≤ _r_ < 1  
 _φ_(_r_) = 0 for 1 ≤ _r_
 
with _r_ a positive distance. Assume that _n_ nearest _from_-points will be used in the interpolation.
An unknown function _f_(**x**) can then be approximated as the weighted sum of _n_ shifted radial basis functions:

_f_(**x**) ≈ Σ<sub>j</sub> _α_<sub>j</sub> _φ_(||**x** − **x**<sub>j</sub>||)

To determine the coefficients _α_<sub>j</sub>, we require that the exact function value is returned at the _n_ _from_-points.
This gives us _n_ equations

_f_(**x**<sub>i</sub>) = _f_<sub>i</sub> 
= Σ<sub>j</sub> _α_<sub>j</sub> _φ_(||**x**<sub>i</sub> − **x**<sub>j</sub>||)

which can be written in matrix form as

**f** = **Φ** · **α**

with **f**, **α** ∈ R<sup>n×1</sup>, and **Φ** ∈ R<sup>n×n</sup>. This system can be solved for the weights-vector **α**.

However, in our case, the _from_-point values vector **f** is not known in advance: it contains the values of the `Variable` that will be interpolated. 

Therefore, the approximation to calculate the interpolatoin in the _to_-point is rewritten as follows:

_f_(**x**<sub>to</sub>) = Σ<sub>j</sub> _α_<sub>j</sub> _φ_(||**x**<sub>to</sub> − **x**<sub>j</sub>||) = **Φ**<sup>T</sup><sub>to</sub> · **α** = **Φ**<sup>T</sup><sub>to</sub> · **Φ**<sup>-1</sup> · **f** = **c**<sup>T</sup> · **f**

The coefficients vector **c** can now be calculated based only on the coordinates by solving the system

**Φ** · **c** = **Φ**<sub>to</sub>.

As every to-point has different nearest neighbours in the _from_-points, the coefficient vector **c** must be calculated for each _to_-point independently. The matrix **Φ** and vector **Φ**<sup>T</sup> must also be calculated for every _to_-point independently.





[//]: # (SOME HANDY MARKDOW URLS:)

[//]: # (MarkDown cheat sheet: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#tables)

[//]: # (render LaTeX eqn as image: https://alexanderrodin.com/github-latex-markdown/)

[//]: # (Greek lower: αβγδεζηϑθικλμνξοπρστυφϕχψω)
[//]: # (Greek upper: ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ)
[//]: # (super, sub: <sup></sup>, <sub></sub> )
[//]: # (operators: + - − · / × √ ∘ ∗)
[//]: # (other: ∂ Δ	∑ ≤ ≥ ∈ )
[//]: # (more HTML math: http://www.unics.uni-hannover.de/nhtcapri/mathematics.html)
[//]: # (and even more: http://www.alanflavell.org.uk/unicode/unidata22.html)