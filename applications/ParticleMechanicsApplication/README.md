# Particle Mechanics Application

This application features continuum-based meshfree and particle methods with main motivations of simulating non-linear large deformable materials, such as free-surface flows, geomechanics phenomena, and extreme events involving impact, penetration, fragmentation, blast, multi-phase interaction, failure evolution, etc.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Readme_files/ParticleMechanicsApplication.gif" width="618" height="280"/>
</p>


The recent research and development have been focused solely on the Material Point Method (MPM. However, the implementation of other continuum-based methods, e.g. the Smoothed Particle Hydrodynamics (SPH), the Galerkin Meshfree Method (GMM), etc, is planned to be done within the same application in the future.

## Getting Started

This application is part of the Kratos Multiphysics Platform. Instructions on how to download, install and run the software in your local machine for development and testing purposes are available for both [Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems.

### Prerequisites

Build [Kratos](https://github.com/KratosMultiphysics/Kratos/wiki) and make sure the following lines are written:

``` cmake
-DEXTERNAL_SOLVER_APPLICATION=ON
-DPARTICLE_MECHANICS_APPLICATION=ON
```

between the compilation configurations, so the Particle Mechanics application is compiled along with auxiliary external solvers required.

## Theory

Particle or meshfree methods are a category of methods where the state of a system is represented by a set of discrete particles, without a fixed connectivity; hence, making such methods suitable for the analysis of moving discontinuities and large deformations with breaking and fragmentation. This approach does not suffer from the mesh distortion and entanglement issues posed by other Lagrangian discretizations such as the finite element method.

### Material Point Method

The MPM is a one of the Lagrangian particle methods which has gained a remarkably increasing popularity due to its capability in simulating solid mechanics problems involving historically dependent materials and large deformations. As MPM is able to combine the strengths of Eulerian and Lagrangian methods, it has been utilized in various engineering applications and industrial purposes, in particular in geomechanics and environmental fluid dynamics field. The method stores the historically changing variables and the material information at the moving particles, the so-called *material points* (MP), and uses a constantly-reset *background mesh* to solve the linear system of equations.

Recommended references for implementation details of MPM in Kratos:
- Iaconeta, I., Larese, A., Rossi, R., & Zhiming, G. (2016). Comparison of a material point method and a Galerkin meshfree method for the simulation of cohesive-frictional materials. *Materials*, 10(10), p. 1150.
- Iaconeta, I., Larese, A., Rossi, R., Oñate, E. (2017). An implicit material point method applied to granular flows. *Procedia Engineering: Proceeding of the 1st International Conference on the Material Point Method (MPM2017)*, Vol 175, 226-232 MPM 2017, Delft, Netherlands.
- Iaconeta, I., Larese, A., Rossi, R., & Oñate, E. (2018). A stabilized mixed implicit Material Point Method for non-linear incompressible solid mechanics. *Computational Mechanics*, 1-18.
- Chandra, B., Larese, A., Iaconeta, I., Rossi, R., Wüchner, R. (2018). Soil-Structure Interaction Simulation of Landslides Impacting a Structure Using an Implicit Material Point Method. *Accepted for publication by Proceeding of the 2nd International Conference on The Material Point Method (MPM2019)*.

## Features

The following features are currently available and subject to development within the Particle Mechanics Application:

- A set of *Neumann* conditions:
    * Static Point load (loads applied directly on the background nodes)
    * Static Line load (a distributed load applied over a line)
    * Static Surface load (a distributed load applied over a face)

- Solid (background) elements:
    * Updated Lagrangian elements - triangular (2D) and tetrahedral (3D), structured and unstructured
    * Updated Lagrangian UP elements - triangular (2D) and tetrahedral (3D), structured and unstructured, with Mixed Variational Methods of displacement and pressure
    * Updated Lagrangian quadrilateral elements - quadrilateral (2D) and hexahedral (3D), structured and unstructured

- Constitutive laws:
    * Linear elastic materials - plane strain, plane stress and 3D
    * Hyperelastic Neo-Hookean laws - finite strain, plane strain and 3D
    * Elasto-plastic laws:
        * Mohr Coulomb - finite strain, associative and non-associative, plane strain and 3D
        * Mohr Coulomb with Strain Softening - finite strain, associative and non-associative, plane strain and 3D
    * Critical state laws:
        * Modified Cam-Clay - finite strain, plane strain and 3D

- Strategies and schemes:
    * Implicit - Newmark prediction and correction scheme

- Other features:
    * Particle erase features - to delete particle outside the interest domain
    * Arbitrary slip boundary condition
    
Some unit tests of the above features can be found in the [tests](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/ParticleMechanicsApplication/tests) folder.

## Available Interfaces

### GiD Interface
It is located in GiD interface repository in [GiD interface repository](https://github.com/KratosMultiphysics/GiDInterface/tree/master/).

Requires [GiD](https://www.gidhome.com/) - Pre and Post Processing software.

## License

The Particle Mechanics Application is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Contact

* **Antonia Larese** - *Group Leader* - [antoldt@cimne.upc.edu](mailto:antoldt@cimne.upc.edu)
* **Bodhinanda Chandra** - *Application Maintainer* - [bodhinanda.chandra@tum.de](mailto:bodhinanda.chandra@tum.de)
