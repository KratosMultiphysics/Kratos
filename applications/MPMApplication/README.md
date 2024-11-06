# MPM Application

- [Overview](#overview)
- [Theory](#theory)
- [Getting Started](#getting-started)
    - [Getting Binaries with `pip` (users)](#getting-binaries-with-pip-users)
    - [Build and Compile Source Code (developers)](#build-and-compile-source-code-developers)
- [Examples](#examples)
- [GiD Interface](#gid-interface)
- [Features](#features)
- [References](#references)
- [License](#license)
- [Contact](#contact)

## Overview

This application implements the **Material Point Method (MPM)** with main motivations of simulating non-linear large deformable materials, such as free-surface flows, geomechanical phenomena, and extreme events involving impact, penetration, fragmentation, blast, multi-phase interaction, failure evolution, etc.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Readme_files/MPMApplication.gif" width="618" height="280"/>
</p>

## Theory

Particle or meshfree methods are a category of methods where the state of a system is represented by a set of particles, without a fixed connectivity; hence, making such methods suitable for the analysis of moving discontinuities and large deformations with breaking and fragmentation. This approach does not suffer from the mesh distortion and entanglement issues posed by other Lagrangian discretizations such as the finite element method.

The **Material Point Method** (MPM) is an hybrid thechnique which uses a fixed background grid (or mesh) for solving the governing equations in a FEM fashion and  set of material particles (MP) for storing all the hystorical variables and material informations. MPM has gained a remarkably increasing popularity due to its capability in simulating  problems involving historically dependent materials and large deformations. As MPM is able to combine the strengths of Eulerian and Lagrangian methods, it has been utilized in various engineering applications and industrial purposes, in particular in geomechanics and environmental fluid dynamics field.

## Getting Started

This application is part of the Kratos Multiphysics framework and it can be obtained either by installing the Kratos binaries with `pip` (suggested for users that want to use the application like a black-box) or by downloading the source code and compiling it (suggested for developers).

### Getting Binaries with `pip` (users)

Kratos binaries are available for Linux, Windows and MacOS and can be obtained with `pip`. Open the terminal and run the following command:

```bash
python3 -m pip install KratosMPMApplication
```

This command will install `KratosMultiphysics` (Kratos Multiphysics Core), `KratosMPMApplication` (application implementing MPM) and `KratosLinearSolversApplication` (dependency of MPMApplication).

### Build and Compile Source Code (developers)

Instructions on how to download, install and run the software in your local machine for development and testing purposes are available for Linux, Windows and MacOS distributions in the [installation page](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md).

In particular, be sure to compile the `MPMApplication` and the auxiliary `LinearSolversApplication`:

* in **Linux**, check that in the `/path_to_kratos/scripts/standard_configure.sh` the following lines are written:
    ```bash
    export KRATOS_APPLICATIONS=
    add_app ${KRATOS_APP_DIR}/MPMApplication
    add_app ${KRATOS_APP_DIR}/LinearSolversApplication
    ```

* in **Windows**, check that in the `/path_to_kratos/scripts/standard_configute.bat` the following lines appears:
    ```console
    set KRATOS_APPLICATIONS=
    CALL :add_app %KRATOS_APP_DIR%\MPMApplication;
    CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
    ```

## Examples
Some use-cases and validation examples are available in the MPM section of the [Examples](https://kratosmultiphysics.github.io/Examples/) repository. Also, some unit tests of the main features can be found in the [tests](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication/tests) folder.

## GiD Interface

A GUI (Graphic User Interface) for the MPM application is also available within the pre and post processing software [GiD](https://www.gidhome.com/). Instructions on how to download and install it are located in the `GiDInterface` [GitHub repository](https://github.com/KratosMultiphysics/GiDInterface/tree/master/). A basic knowledge of GiD is required.

## Features

The following features are currently available and subject to development within the `MPMApplication`.

**Formulations**
* Irreducible formulation (u displacement based)
* Mixed UP (displacement/pressure) formulation

**Element types**
* Updated Lagrangian elements - triangular and quadrilateral (2D) and tetrahedral and hexahedral (3D), structured and unstructured, using classical or partitioned quadrature rules (this latter limited to explicit MPM)
* Updated Lagrangian axis-symmetric elements - triangular and quadrilateral (2D), structured and unstructured
* Updated Lagrangian mixed UP elements - triangular (2D) and tetrahedral (3D), structured and unstructured, stabilized using  Variational Multiscale Stabilization (VMS) or Pressure Projection techniques

**Constitutive laws**
* Linear isotropic elastic materials - plane strain, plane stress, axis-symmetric, and 3D
* Hyperelastic Neo-Hookean laws - finite strain, plane strain, axis-symmetric, and 3D
* Elasto-plastic laws:
    * Mohr Coulomb - finite strain, associative and non-associative, plane strain, axis-symmetric, and 3D
    * Mohr Coulomb with Strain Softening - finite strain, associative and non-associative, plane strain, axis-symmetric, and 3D
* Critical state laws:
    * Modified Cam-Clay - finite strain, plane strain, axis-symmetric, and 3D
    * Johnson Cook Thermal Plastic (just for explicit MPM)

**Boundary conditions**
* Grid-Based Conditions (conforming): applied directly at the background nodes
    * Neumann: Point load
    * Neumann: Line load (a distributed load applied over a line)
    * Neumann: Surface load (a distributed load applied over a face)
    * Dirichlet: Slip and non-slip conditions for arbitrary inclination
* Material Point-Based Conditions (non-conforming): applied on movable boundary particles
    * Neumann:
        * moving point load
        * interface condition for partitioned coupling with DEM
    * Dirichlet: fixed, slip or contact condition
        * penalty method
        * Lagrange multiplier method
        * perturbed Lagrangian method
        * interface condition for partitioned coupling with FEM, RBS,...

**Time schemes**
* Implicit - Newmark/Bossak prediction and correction scheme for static, quasi-static, and dynamic problems
* Explicit

**Other features**
* Partitioned coupling with Finite Element Method (FEM) - weak and strong coupling of nonconforming discretization
* Partitioned coupling with the Discrete Element Method (DEM)
* Partitioned coupling with the Rigid Body Solver (RBS)
* Material point erase features - to delete material points outside the interest domain


## References

Recommended references for implementation details of MPM in Kratos:

* Singer, V., (2024). **Partitioned Coupling Strategies to Simulate the Impact of Granular Mass Flows on Flexible Protective Structures**, *PhD Thesis*, Technical University of Munich.
* Singer, V., Teschemacher, T., Larese, A., Wüchner, R., Bletzinger, K.U. (2024). **Lagrange multiplier imposition of non-conforming essential boundary conditions in implicit Material Point Method**, *Computational Mechanics*, 73, 1311–1333 DOI: <a href="https://doi.org/10.1007/s00466-023-02412-w">10.1007/s00466-023-02412-w</a>.
* Singer, V., Sautter, K.B., Larese, A., Wüchner, R., Bletzinger K.-U., (2023) **Partitioned Coupling Approaches for the Simulation of Natural Hazards Impacting Protective Structures**, *VIII International Conference on Particle-Based Methods*. DOI: <a href="https://doi.org/10.23967/c.particles.2023.002">10.23967/c.particles.2023.002</a>.
* Singer, V., Larese, A., Wüchner, R., Bletzinger K.-U., (2023). **Partitioned MPM-FEM Coupling Approach for Advanced Numerical Simulation of Mass-Movement Hazards Impacting Flexible Protective Structures**, *X International Conference on Computational Methods for Coupled Problems in Science and Engineering*. DOI: <a href="https://doi.org/10.23967/c.coupled.2023.026">10.23967/c.coupled.2023.026</a>.
* Singer, V., Sautter, K.B., Larese, A., Wüchner, R., Bletzinger, K.-U. (2022). **A partitioned material point method and discrete element method coupling scheme**, *Advanced Modeling and Simulation in Engineering Sciences*, 9(16). DOI: <a href="https://doi.org/10.1186/s40323-022-00229-5">doi.org/10.1186/s40323-022-00229-5</a>.
* Wilson, P., (2022). **A computational impact analysis approach leveraging non-conforming spatial, temporal and methodological discretisations**, *PhD Thesis*, University of Queensland. DOI: <a href="https://doi.org/10.14264/3e10f66">10.14264/3e10f66</a>.
* Singer, V., Bodhinanda, C., Larese, A., Wüchner, R., Bletzinger K.-U., (2021). **A Staggered Material Point Method and Finite Element Method Coupling Scheme Using Gauss Seidel Communication Pattern**, *9th edition of the International Conference on Computational Methods for Coupled Problems in Science and Engineering*. DOI: <a href="https://doi.org/10.23967/coupled.2021.006">10.23967/coupled.2021.006</a>.
* Chandra, B., Singer, V., Teschemacher, T., Wuechner, R., & Larese, A. (2021). **Nonconforming Dirichlet boundary conditions in implicit material point method by means of penalty augmentation**, *Acta Geotechnica*, 16(8), 2315-2335. DOI: <a href="https://doi.org/10.1007/s11440-020-01123-3">10.1007/s11440-020-01123-3</a>.
* Wilson, P., Wüchner, R., & Fernando, D. (2021). **Distillation of the material point method cell crossing error leading to a novel quadrature‐based C0 remedy**, *International Journal for Numerical Methods in Engineering*, 122(6), 1513-1537. DOI: <a href="https://doi.org/10.1002/nme.6588">10.1002/nme.6588</a>.
* Iaconeta, I., Larese, A., Rossi, R., & Oñate, E. (2018). **A stabilized mixed implicit Material Point Method for non-linear incompressible solid mechanics**, *Computational Mechanics*, 1-18. DOI: <a href="https://doi.org/10.1007/s00466-018-1647-9">10.1007/s00466-018-1647-9</a>.
* Iaconeta, I., Larese, A., Rossi, R., & Zhiming, G. (2016). **Comparison of a material point method and a Galerkin meshfree method for the simulation of cohesive-frictional materials**, *Materials*, 10(10), p. 1150. DOI: <a href="https://doi.org/10.3390/ma10101150">10.3390/ma10101150</a>.

## License

The `MPMApplication` is **OPEN SOURCE**. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The **BSD** (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Contact

* **Antonia Larese** - *Group Leader* - [antonia.larese@unipd.it](mailto:antonia.larese@unipd.it)
* **Veronika Singer** - *Developer* - [veronika.singer@tum.de](mailto:veronika.singer@tum.de)
* **Laura Moreno** - *Developer* - [laura.morenomartinez@ua.es](mailto:laura.morenomartinez@ua.es)
* **Andi Makarim Katili** - *Developer* - [andi.katili@tum.de](mailto:andi.katili@tum.de)
