# MPM Application

This application implements the Material Point Method (MPM) with main motivations of simulating non-linear large deformable materials, such as free-surface flows, geomechanical phenomena, and extreme events involving impact, penetration, fragmentation, blast, multi-phase interaction, failure evolution, etc.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Readme_files/MPMApplication.gif" width="618" height="280"/>
</p>


## Getting Started

This application is part of the Kratos Multiphysics Platform. Instructions on how to download, install and run the software in your local machine for development and testing purposes are available for both Linux and Windows distributions [Installation page](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md).

### Prerequisites

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

### GiD Interface
A GiD user interface for the MPM application is also available. It is located in GiD interface repository in [GiD interface repository](https://github.com/KratosMultiphysics/GiDInterface/tree/master/).

It requires [GiD](https://www.gidhome.com/) - Pre and Post Processing software.

## Theory

Particle or meshfree methods are a category of methods where the state of a system is represented by a set of particles, without a fixed connectivity; hence, making such methods suitable for the analysis of moving discontinuities and large deformations with breaking and fragmentation. This approach does not suffer from the mesh distortion and entanglement issues posed by other Lagrangian discretizations such as the finite element method.

### Material Point Method

The MPM is an hybrid thechnique which uses a fixed background grid (or mesh) for solving the governing equations in a FEM fashion and  set of material particles (MP) for storing all the hystorical variables and material informations. MPM has gained a remarkably increasing popularity due to its capability in simulating  problems involving historically dependent materials and large deformations. As MPM is able to combine the strengths of Eulerian and Lagrangian methods, it has been utilized in various engineering applications and industrial purposes, in particular in geomechanics and environmental fluid dynamics field.

Recommended references for implementation details of MPM in Kratos:
- Singer, V.; Partitioned Coupling Strategies to Simulate the Impact of Granular Mass Flows on Flexible Protective Structures; PhD Thesis, Technical University of Munich (2024)
- Singer, V.; Teschemacher T.; Larese A.; Wüchner R. Bletzinger K.-U.; Lagrange multiplier imposition of non-conforming essential boundary conditions in implicit material point method. Comput Mech 73, 1311–1333 (2024). https://doi.org/10.1007/s00466-023-02412-w
- Singer, V.; Sautter, K.B.; Larese, A.; Wüchner, R.; Bletzinger K.-U.; Partitioned Coupling Approaches for the Simulation of Natural Hazards Impacting Protective Structures. VIII International Conference on Particle-Based Methods (2023); doi: https://doi.org/10.23967/c.particles.2023.002
- Singer, V.; Larese, A.; Wüchner, R.; Bletzinger K.-U.; Partitioned MPM-FEM Coupling Approach for Advanced Numerical Simulation of Mass-Movement Hazards Impacting Flexible Protective Structures,  X International Conference on Computational Methods for Coupled Problems in Science and Engineering (2023); doi: https://doi.org/10.23967/c.coupled.2023.026
- Singer, V.; Sautter, K.B., Larese, A., Wüchner, R.; Bletzinger, K.U.; A partitioned material point method and discrete element method coupling scheme, Advanced Modeling and Simulation in Engineering Sciences, 9(16), (2022); DOI: https://doi.org/10.1186/s40323-022-00229-5
- Wilson, P.; A computational impact analysis approach leveraging non-conforming spatial, temporal and methodological discretisations; PhD Thesis, University of Queensland (2022); doi: https://doi.org/10.14264/3e10f66
- Singer, V.; Bodhinanda, C.; Larese, A.; Wüchner, R.; Bletzinger K.-U.; A Staggered Material Point Method and Finite Element Method Coupling Scheme Using Gauss Seidel Communication Pattern, 9th edition of the International Conference on Computational Methods for Coupled Problems in Science and Engineering (2021); doi: https://doi.org/10.23967/coupled.2021.006
- Chandra, B., Singer, V., Teschemacher, T., Wuechner, R., & Larese, A. (2021). Nonconforming Dirichlet boundary conditions in implicit material point method by means of penalty augmentation. Acta Geotechnica, 16(8), 2315-2335. DOI: https://doi.org/10.1007/s11440-020-01123-3
- Wilson, P., Wüchner, R., & Fernando, D. (2021). Distillation of the material point method cell crossing error leading to a novel quadrature‐based C 0 remedy. International Journal for Numerical Methods in Engineering, 122(6), 1513-1537, doi:https://doi.org/10.1002/nme.6588
- Iaconeta, I., Larese, A., Rossi, R., & Oñate, E. (2018). A stabilized mixed implicit Material Point Method for non-linear incompressible solid mechanics. *Computational Mechanics*, 1-18. DOI https://doi.org/10.1007/s00466-018-1647-9
- Iaconeta, I., Larese, A., Rossi, R., & Zhiming, G. (2016). Comparison of a material point method and a Galerkin meshfree method for the simulation of cohesive-frictional materials. *Materials*, 10(10), p. 1150. doi: https://doi.org/10.3390/ma10101150
-
## Features

The following features are currently available and subject to development within the MPM Application:
- Formulation:
  * Irreducible formulations (u displacement based)
  * Mixed UP formulations

- Element types:
    * Updated Lagrangian elements - triangular and quadrilateral (2D) and tetrahedral and hexahedral (3D), structured and unstructured, using classical or partitioned quadrature rules (this latter limited to explicit MPM)
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
    * Grid-Based Conditions (conforming): applied directly at the background nodes
        * Neumann: Point load
        * Neumann: Line load (a distributed load applied over a line)
        * Neumann: Surface load (a distributed load applied over a face)
        * Dirichlet: Slip and non-slip conditions for arbitrary inclination.
    * Material Point-Based Conditions (non-conforming): applied on movable boundary particles
        * Neumann: 
            * moving point load 
            * interface condition for partitioned coupling with DEM
        * Dirichlet: fixed, slip or contact condition
            * penalty method, 
            * Lagrange multiplier method 
            * perturbed Lagrangian method
            * interface condition for partitioned coupling with FEM, RBS,...

- Strategies and schemes:
    * Implicit - Newmark/Bossak prediction and correction scheme for static, quasi-static, and dynamic problems
    * Explicit

- Other features:
    * Partitioned coupling with Finite Element Method (FEM) - weak and strong coupling of nonconforming discretization
    * Partitioned coupling with the Discrete Element Method (DEM)
    * Partitioned coupling with the Rigid Body Solver (RBS)
    * material point erase features - to delete material points outside the interest domain

## License

The MPM Application is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Contact

* **Antonia Larese** - *Group Leader* - [antonia.larese@unipd.it](mailto:antonia.larese@unipd.it)
* **Veronika Singer** - *Developer* - [veronika.singer@tum.de](mailto:veronika.singer@tum.de)
* **Laura Moreno** - *Developer* - [laura.morenomartinez@unipd.it](mailto:laura.morenomartinez@unipd.it)
* **Andi Makarim Katili** - *Developer* - [andi.katili@tum.de](mailto:andi.katili@tum.de)
