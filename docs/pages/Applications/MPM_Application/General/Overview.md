---
title: MPM Application
keywords: 
tags: [Overview.md]
sidebar: mpm_application
summary: 
---

This application implements the Material Point Method (MPM) with main motivations of simulating non-linear large deformable materials, such as free-surface flows, geomechanical phenomena, and extreme events involving impact, penetration, fragmentation, blast, multi-phase interaction, failure evolution, etc.

![MPMApplication](https://user-images.githubusercontent.com/51473791/191960884-1f1c5a0c-efec-40ca-ac6d-2d53b5530739.gif)

## Getting Started

This application is part of the Kratos Multiphysics Platform. Instructions on how to download, install and run the software in your local machine for development and testing purposes are available for both Linux and Windows distributions [Installation page](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md).

**Prerequisites**

Build Kratos and check the [configuration files](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md#configuration-scripts-examples)

In LINUX: check that in the /path_to_kratos/scripts/configure.sh the followinglines are written:

``` cmake
-DMPM_APPLICATION=ON
-DLINEAR_SOLVERS_APPLICATION=ON
```

In WINDOWS: check that in the /path_to_kratos/scripts/configute.bat the following lines appears: 

```set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\MPMApplication;
```

so the MPM application is compiled along with auxiliary linear solvers required.

## Examples
Some use-cases and validation examples are available in the MPM section of the [Examples](https://kratosmultiphysics.github.io/Examples/) repository. Also, some unit tests of the main features can be found in the [tests](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication/tests) folder.

**GiD Interface**

A GiD user interface for the MPM application is also available. It is located in GiD interface repository in [GiD interface repository](https://github.com/KratosMultiphysics/GiDInterface/tree/master/).
It requires [GiD](https://www.gidhome.com/) - Pre and Post Processing software.

## Theory

Particle or meshfree methods are a category of methods where the state of a system is represented by a set of particles, without a fixed connectivity; hence, making such methods suitable for the analysis of moving discontinuities and large deformations with breaking and fragmentation. This approach does not suffer from the mesh distortion and entanglement issues posed by other Lagrangian discretizations such as the finite element method.

**Material Point Method**

The MPM is an hybrid thechnique which uses a fixed background grid (or mesh) for solving the governing equations in a FEM fashion and  set of material particles (MP) for storing all the hystorical variables and material informations. MPM has gained a remarkably increasing popularity due to its capability in simulating  problems involving historically dependent materials and large deformations. As MPM is able to combine the strengths of Eulerian and Lagrangian methods, it has been utilized in various engineering applications and industrial purposes, in particular in geomechanics and environmental fluid dynamics field.

Recommended references for implementation details of MPM in Kratos:
- Singer, V.; Sautter, K.B., Larese, A., Wüchner, R.; Bletzinger, K.U.; A partitioned material point method and discrete element method coupling scheme, Advanced Modeling and Simulation in Engineering Sciences, 9(16), (2022); DOI: https://doi.org/10.1186/s40323-022-00229-5
- Chandra, B., Singer, V., Teschemacher, T., Wuechner, R., & Larese, A. (2021). Nonconforming Dirichlet boundary conditions in implicit material point method by means of penalty augmentation. Acta Geotechnica, 16(8), 2315-2335. DOI: 10.1007/s11440-020-01123-3
- Wilson, P., Wüchner, R., & Fernando, D. (2021). Distillation of the material point method cell crossing error leading to a novel quadrature‐based C 0 remedy. International Journal for Numerical Methods in Engineering, 122(6), 1513-1537.
- Iaconeta, I., Larese, A., Rossi, R., & Oñate, E. (2018). A stabilized mixed implicit Material Point Method for non-linear incompressible solid mechanics. *Computational Mechanics*, 1-18. DOI 10.1007/s00466-018-1647-9
- Iaconeta, I., Larese, A., Rossi, R., & Zhiming, G. (2016). Comparison of a material point method and a Galerkin meshfree method for the simulation of cohesive-frictional materials. *Materials*, 10(10), p. 1150. doi: 10.3390/ma10101150


## Features

The following features are currently available and subject to development within the MPM Application:
- Formulation:
  * Irreducible formulations (U displacement based)
  * Mixed UP formulations

- Element types:
    * Updated Lagrangian elements - triangular and quadrilatera (2D) and tetrahedral and hexahedral (3D), structured and unstructured, using classical or partitioned quadrature rules (this latter limited to explicit MPM)
    * Updated Lagrangian axis-symmetric elements - triangular and quadrilateral (2D), structured and unstructured
    * Updated Lagrangian mixed UP elements - triangular (2D) and tetrahedral (3D), structured and unstructured, stabilized using  Variational Multiscale Stabilization or Pressure Projection techniques

- Constitutive laws:
    * Linear isotropic elastic materials - plane strain, plane stress, axis-symmetric, and 3D
    * Hyperelastic Neo-Hookean laws - finite strain, plane strain, axis-symmetric, and 3D
    * Elasto-plastic laws:
        * Mohr Coulomb - finite strain, associative and non-associative, plane strain, axis-symmetric, and 3D
        * Mohr Coulomb with Strain Softening - finite strain, associative and non-associative, plane strain, axis-symmetric, and 3D
    * Critical state laws:
        * Modified Cam-Clay - finite strain, plane strain, axis-symmetric, and 3D
        * Johnson Cook Thermal Plastic (just for explicit MPM)

- A set of Boundary conditions:
    * Grid-Based Conditions (conforming): applied directly on the background nodes
        * Neumann: Point load
        * Neumann: Line load (a distributed load applied over a line)
        * Neumann: Surface load (a distributed load applied over a face)
        * Dirichlet: Slip and No slip condition on arbitrary boundary.
    * Material Point-Based Conditions (non conforming): applied in moveable boundary particles
        * Neumann: Moving point load
        * Dirichlet: Imposition of displacements (homogeneous and inhomogeneous) using penalty method

- Strategies and schemes:
    * Implicit - Newmark/Bossak prediction and correction scheme for static, quasi-static, and dynamic problems
    * Explicit

- Other features:
    * Partitioned coupling with Finite Element Method - weak and strong coupling of nonconforming discretization
    * Partitioned coupling with the Discrete Element Method
    * material point erase features - to delete material points outside the interest domain


## License

The MPM Application is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Contact

* **Antonia Larese** - *Group Leader* - [antonia.larese@unipd.it](mailto:antonia.larese@unipd.it)
* **Veronika Singer** - *Developer* - [veronika.singer@tum.de](mailto:veronika.singer@tum.de)
* **Laura Moreno** - *Developer* - [laura.morenomartinez@unipd.it](mailto:laura.morenomartinez@unipd.it)
