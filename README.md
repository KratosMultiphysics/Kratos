<img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/kratos.png">

[![Release][release-image]][releases] [![License][license-image]][license] [![Master][kratos-master-status]][travis-branches]

_KRATOS Multiphysics_ ("Kratos") is a framework for building parallel, multi-disciplinary simulation software, aiming at modularity, extensibility, and high performance. Kratos is written in C++, and counts with an extensive Python interface. More in [Overview](https://github.com/KratosMultiphysics/Kratos/wiki/Overview)

**Kratos** is **free** under BSD-4 [license](https://github.com/KratosMultiphysics/Kratos/wiki/Licence) and can be used even in comercial softwares as it is. Many of its main applications are also free and BSD-4 licensed but each derived application can have its own propietary license. More in [application licensing model](link to application licensising model webpage)

[release-image]: https://img.shields.io/badge/release-5.2-green.svg?style=flat
[releases]: https://github.com/KratosMultiphysics/Kratos/releases

[license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/license.txt
[kratos-master-status]: https://travis-ci.org/KratosMultiphysics/Kratos.svg?branch=master
[travis-branches]: https://travis-ci.org/KratosMultiphysics/Kratos/branches

# Main Features
**Kratos** is multiplatform and available for Windows,Linux (several distros) and can be compiled in OSX.

**Kratos** is OpenMP and MPI parallel and scalable up to thousands of cores.

**Kratos** provides a core which defines the common framework and several application which work like plug-ins and extend it in different fields.

Its main applications are:
- [DEM]() for cohesive and non cohesive shperic and non spheric particles simultion
- [Fluid Dynamics]() Provides 2D and 3D incompressible fluids formulation
- [Fluid Structure Interaction]() for solution of different FSI problems
- [Structural Mechanics]() Providing solution for solid, shell and beam structures with linear and nonlinear, static and dynamic behavior

Some main modules are:
- [External Solvers]()
- [Trilinos]()
- [Metis]()
- [Meshing]()

# Examples of use
Kratos has been used for simulation of many different problems in a wid variaty of disciplines ranging from wind over singular building to granular domain dynamics. Some examples and validation benchmarks simulated by Kratos can be found [here](https://kratosmultiphysics.github.io/Examples/)

# Contributors
Organizations contributing to Kratos: 
<table>
<tr>
  <td><img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Logos/CIMNE_logo.png" width="128"></td>
  <td><img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Logos/TUM_Logo.png" width="128"></td>
</tr>
<tr>
  <td>International Center for Numerical Methods in Engineering</td>
  <td>Chair of Structural Analysis<br>
Technical University of Munich
</td>
</tr>
</table>

  
# Known Users
Some users of the technologies developed in Kratos are:

<table>
<tr>
  <td><img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Logos/AIRBUS_logo.png" width="128"></td>
  <td><img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Logos/siemens_logo.png" width="128"></td>
  <td><img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Logos/onera_logo.png" width="128"></td>
</tr>
<tr>
  <td>Airbus Defence and Space<br>Stress Methods & Optimisation Department</td>
  <td>Siemens AG<br>
Corporate Technology
</td>
  <td>ONERA, The French Aerospace Lab<br>
Applied Aerodynamics Department 

</td>
</tr>
</table>

Looking forward to seeing your logo here! 

# Special Thanks To
In Kratos Core:
- [Boost](http://www.boost.org/) for boost.python and ublas
- [GidPost](https://www.gidhome.com/gid-plus/tools/476/gidpost/) providing output to [GiD](https://www.gidhome.com/)
- [AMGCL](https://github.com/ddemidov/amgcl) for its highly scalable multigrid solver
- ZLib 

In applications
- [Trilinos](https://trilinos.org/) for MPI linear algebra and solvers used in trilinos application
- [METIS](http://glaros.dtc.umn.edu/gkhome/views/metis) for partitioning in metis application
