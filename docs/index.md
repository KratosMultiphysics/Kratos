---
keywords: Kratos Main Overview Readme
tags: []
sidebar: kratos_sidebar
permalink: index.html
---

<p align=center><img height="72.125%" width="72.125%" src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/kratos.png"></p>

[![Release][release-image]][releases] [![License][license-image]][license] [![Github CI][Nightly-Build]][Nightly-link] [![DOI][DOI-image]][DOI]

[release-image]: https://img.shields.io/badge/release-8.1-green.svg?style=flat
[releases]: https://github.com/KratosMultiphysics/Kratos/releases

[license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/license.txt

[Nightly-Build]: https://github.com/KratosMultiphysics/Kratos/workflows/Nightly%20Build/badge.svg
[Nightly-link]: https://github.com/KratosMultiphysics/Kratos/actions?query=workflow%3A%22Nightly+Build%22

[DOI-image]: https://zenodo.org/badge/DOI/10.5281/zenodo.3234644.svg

[DOI]: https://doi.org/10.5281/zenodo.3234644

_KRATOS Multiphysics_ ("Kratos") is a framework for building parallel, multi-disciplinary simulation software, aiming at modularity, extensibility, and high performance. Kratos is written in C++, and counts with an extensive Python interface. More in [Overview](pages/Overview)

**Kratos** is **free** under BSD-4 [license](pages/Licence) and can be used even in comercial softwares as it is. Many of its main applications are also free and BSD-4 licensed but each derived application can have its own propietary license.


# Main Features
**Kratos** is __multiplatform__ and available for __Windows, Linux__ (several distros) and __macOS__.

**Kratos** is __OpenMP__ and __MPI__ parallel and scalable up to thousands of cores.

**Kratos** provides a core which defines the common framework and several application which work like plug-ins that can be extended in diverse fields.

Its main applications are:
- [DEM](applications/DEMApplication) for cohesive and non cohesive spheric and non spheric particles simultion
- [Fluid Dynamics](applications/FluidDynamicsApplication/README.md) Provides 2D and 3D incompressible fluids formulation
- [Fluid Structure Interaction](applications/FSIApplication/README.md) for solution of different FSI problems
- [Structural Mechanics](applications/StructuralMechanicsApplication/README.md) Providing solution for solid, shell and beam structures with linear and nonlinear, static and dynamic behavior
- [Contact Structural Mechanics](applications/ContactStructuralMechanicsApplication/README.md) For contact problems used along the [Structural Mechanics application](applications/StructuralMechanicsApplication/README.md)

Some main modules are:
- [Linear Solvers](applications/LinearSolversApplication/README.md)
- [Trilinos](applications/TrilinosApplication/README.md)
- [Metis](applications/MetisApplication/README.md)
- [Meshing](applications/MeshingApplication/README.md)

<!-- We are alreayd on the doc page! -->
<!-- # Documentation
Here you can find the basic documentation of the project:
## Getting Started
* Getting Kratos (Last compiled Release)
    * [Kratos for Linux](Getting-Kratos-Binaries-for-Linux)
    * [Kratos for Windows](Getting-Kratos-Binaries-for-Windows)
    // * [Kratos for Mac](MacOS-install) 
    * [Kratos for GiD](Getting-Kratos-binaries-(via-GiD))
* Compiling Kratos
    * [See INSTALL.md](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md)

## Tutorials
* [Running an example from GiD](Running-an-example-from-GiD)
* [Kratos input files and I/O](Kratos-input-files-and-IO)
* [Data management](Data-management)
* [Solving strategies](Solving-strategies)
* [Manipulating solution values](Manipulating-solution-values)
* [Multiphysics](Multiphysics-example)

## More documentation
[Wiki](https://github.com/KratosMultiphysics/Kratos/wiki) -->

# Examples of use
Kratos has been used for simulation of many different problems in a wide variety of disciplines ranging from wind over singular building to granular domain dynamics. Some examples and validation benchmarks simulated by Kratos can be found [here](https://kratosmultiphysics.github.io/Examples/)

<div class="row">
    <div class="col-md-12" style="text-align:center;">
        <img src="https://github.com/KratosMultiphysics/Examples/raw/master/fluid_dynamics/use_cases/barcelona_wind/resources/BarcelonaVelocityVector.png" width="75%">
        <p>Barcelona Wind Simulation</p>
    </div>
</div>

# Contributors
Organizations contributing to Kratos:

<div class="row">
    <div class="col-md-6" style="text-align:center;">
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/CIMNE_logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>International Center</span>
        <br><span>for Numerical Methods in Engineering</span><br>
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/TUM_Logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>Chair of Structural Analysis</span>
        <br><span>Technical University of Munich</span><br>
    </div>
    <div class="col-md-6" style="text-align:center;">
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/altair-sponsor-logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>Altair Engineering</span><br>
        <br>
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/Deltares_logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>Deltares</span><br>
        <br>
    </div>
</div>

# Our Users
Some users of the technologies developed in Kratos are:

<div class="row">
    <div class="col-md-6" style="text-align:center;">
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/AIRBUS_logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>Airbus Defence and Space</span>
        <br><span>Stress Methods & Optimisation Department</span><br>
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/siemens_logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>Siemens AG</span>
        <br><span>Corporate Technology</span><br>
    </div>
    <div class="col-md-6" style="text-align:center;">
        <img src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/onera_logo.png" style="margin-bottom:15px;max-height:80px;min-height:80px">
        <br><span>ONERA, The French Aerospace Lab</span>
        <br><span>Applied Aerodynamics Department</span><br>
    </div>
</div>

# Special Thanks To
In Kratos Core:
- [Boost](http://www.boost.org/) for ublas
- [pybind11](https://github.com/pybind/pybind11) for exposing C++ to python
- [GidPost](https://www.gidhome.com/gid-plus/tools/476/gidpost/) providing output to [GiD](https://www.gidhome.com/)
- [AMGCL](https://github.com/ddemidov/amgcl) for its highly scalable multigrid solver
- [JSON](https://github.com/nlohmann/json) JSON for Modern C++
- [filesystem](https://github.com/gulrak/filesystem) Header-only single-file std::filesystem compatible helper library, based on the C++17 specs
- [ZLib](https://zlib.net/) The compression library

In applications:
- [Eigen](http://eigen.tuxfamily.org) For linear solvers used in the [LinearSolversApplication](applications/LinearSolversApplication)
- [Trilinos](https://trilinos.org/) for MPI linear algebra and solvers used in [TrilinosApplication](applications/TrilinosApplication)
- [METIS](http://glaros.dtc.umn.edu/gkhome/views/metis) for partitioning in [MetisApplication](applications/MetisApplication/README.md)

# How to cite Kratos?
Please, use the following references when citing Kratos in your work.
- Dadvand, P., Rossi, R. & Oñate, E. An Object-oriented Environment for Developing Finite Element Codes for Multi-disciplinary Applications. Arch Computat Methods Eng 17, 253–297 (2010). https://doi.org/10.1007/s11831-010-9045-2
- Dadvand, P., Rossi, R., Gil, M., Martorell, X., Cotela, J., Juanpere, E., Idelsohn, S., Oñate, E. (2013). Migration of a generic multi-physics framework to HPC environments. Computers & Fluids. 80. 301–309. 10.1016/j.compfluid.2012.02.004.
- Mataix Ferrándiz, V., Bucher, P., Rossi, R., Cotela, J., Carbonell, J. M., Zorrilla, R., … Tosi, R. (2020, November 27). KratosMultiphysics (Version 8.1). Zenodo. https://doi.org/10.5281/zenodo.3234644
