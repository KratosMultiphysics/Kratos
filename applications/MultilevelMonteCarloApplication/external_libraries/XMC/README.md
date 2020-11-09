[DOI-image]: https://zenodo.org/badge/DOI/10.5281/zenodo.3235832.svg

[DOI]: https://doi.org/10.5281/zenodo.3235832

# XMC

[![DOI][DOI-image]][DOI]

**XMC** is a Python library for parallel, adaptive, hierarchical Monte Carlo algorithms, aiming at reliability, modularity, extensibility and high performance.

**XMC** is developed within the [ExaQUte](http://exaqute.eu/) European H2020 project.

**XMC** is **free** under BSD-4 license.

## Main Features
The algorithms **XMC** can run include:
- Monte Carlo,
- Multilevel Monte Carlo,
- Continuation Multilevel Monte Carlo,
- Asynchronous Monte Carlo and Multilevel Monte Carlo.

## Documentation and Usage
Documentation can be found in the form of Docstrings in the code and [here](http://exaqute.eu/wp-content/uploads/sites/10/2019/10/M12_ExaQUte_deliverable_5.2_Release_of_ExaQUte-MLMC_Python_engine.pdf).
Some examples and validation benchmarks can be found [here](https://gitlab.com/RiccardoRossi/exaqute-xmc/-/tree/development/examples/).

## Dependencies
- NumPy;
- SciPy;
- [COMPSs](https://github.com/bsc-wdc/compss) (including its Python interface PyCOMPSs) or [HyperLoom](https://github.com/It4innovations/HyperLoom) for parallel computation (optional).

## External collaborations
**XMC** is integrated with [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos) as solver software.

## How to cite XMC?
All the necessary metadata, as well as formatted citations, are provided on the [Zenodo record](http://doi.org/10.5281/zenodo.3235832).