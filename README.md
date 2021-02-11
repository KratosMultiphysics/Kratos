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

_KRATOS Multiphysics_ ("Kratos") is a framework for building parallel, multi-disciplinary simulation software, aiming at modularity, extensibility, and high performance. Kratos is written in C++, and counts with an extensive Python interface. More in [Overview](https://github.com/KratosMultiphysics/Kratos/wiki/Overview)

**Kratos** is **free** under BSD-4 [license](https://github.com/KratosMultiphysics/Kratos/wiki/Licence) and can be used even in comercial softwares as it is. Many of its main applications are also free and BSD-4 licensed but each derived application can have its own propietary license.


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

# Examples of use
Kratos has been used for simulation of many different problems in a wide variety of disciplines ranging from wind over singular building to granular domain dynamics. Some examples and validation benchmarks simulated by Kratos can be found [here](https://kratosmultiphysics.github.io/Examples/)

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Examples/raw/master/fluid_dynamics/use_cases/barcelona_wind/resources/BarcelonaVelocityVector.png" width="288">
  Barcelona Wind Simulation
</span>
<br>

# Contributors
Organizations contributing to Kratos:

<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/CIMNE_logo.png" width="128">
</br></br><p>International Center for Numerical Methods in Engineering</p>

</br></br>

<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/TUM_Logo.png" width="128">
</br><p>Chair of Structural Analysis</br>Technical University of Munich </p>

<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/altair-sponsor-logo.png" width="128">
</br><p>Altair Engineering</p>

<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/Deltares_logo.png" width="128">
</br><p>Deltares</p>

# Our Users
Some users of the technologies developed in Kratos are:

<span>
<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/AIRBUS_logo.png" width="128">
<p>Airbus Defence and Space</br>Stress Methods & Optimisation Department</p>
</span>
<span>
<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/siemens_logo.png" width="128">
<p>Siemens AG</br>Corporate Technology</p>
</span>
<span>
<img align="left" src="https://github.com/KratosMultiphysics/Documentation/raw/master/Wiki_files/Logos/onera_logo.png" width="128">
<p>ONERA, The French Aerospace Lab<br>Applied Aerodynamics Department</p>
</span>


Looking forward to seeing your logo here!

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
