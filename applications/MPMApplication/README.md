# MPM Application

This application implements the **Material Point Method (MPM)** with main motivations of simulating non-linear large deformable materials, such as free-surface flows, geomechanical phenomena, and extreme events involving impact, penetration, fragmentation, blast, multi-phase interaction, failure evolution, etc.

![MPMApplication](https://user-images.githubusercontent.com/51473791/191960884-1f1c5a0c-efec-40ca-ac6d-2d53b5530739.gif)

## Theory

Particle or meshfree methods are a family of methods in which the state of a system is represented by a set of particles, without a fixed connectivity. As a consequence, these methods are particularly well suited for the analysis of moving discontinuities and large deformations with breaking and fragmentation. This approach does not suffer from the mesh distortion and entanglement issues posed by other Lagrangian discretizations such as the Finite Element Method (FEM).

The **Material Point Method** (MPM) is an hybrid thechnique which uses a fixed background grid (or mesh) for solving the governing equations in a FEM fashion and  set of material particles (MP) for storing all the hystorical variables and material informations. The MPM has gained a remarkably increasing popularity due to its capability in simulating problems involving historically dependent materials and large deformations. As MPM is able to combine the strengths of both Eulerian and Lagrangian methods, it has been used in various engineering applications and industrial purposes, in particular in geomechanics and in the environmental fluid dynamics field.

## Getting Started

The `MPMApplication` is part of the Kratos Multiphysics framework and can be obtained in two different ways:
* by installing the Kratos binaries using the package manager `pip` (suggested for users that want to use the application like a black-box);
* by downloading the source code and compiling it (suggested for developers).

### Getting Binaries with `pip` (users)

Kratos binaries are available for Linux, Windows and MacOS and can be installed by using the `pip` package manager.

Open the terminal and run the following command:

```bash
python3 -m pip install KratosMPMApplication
```

This command will install the following packages:
* `KratosMultiphysics`: Kratos Multiphysics Core;
* `KratosMPMApplication`: application implementing MPM;
* `KratosLinearSolversApplication`: dependency required by `MPMApplication`.

### Build from Source (developers)

Instructions on how to download, compile and run Kratos in your local machine for development and testing purposes are available for Linux, Windows and MacOS distributions in the [installation page](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md).

In particular, in order to use the `MPMApplication` it is also required to compile the auxiliary `LinearSolversApplication`.

* In **Linux**, the following lines must appear in the `/path_to_kratos/scripts/standard_configure.sh` file:
    ```bash
    export KRATOS_APPLICATIONS=
    add_app ${KRATOS_APP_DIR}/MPMApplication
    add_app ${KRATOS_APP_DIR}/LinearSolversApplication
    ```

* In **Windows**, the following lines must appear in the `/path_to_kratos/scripts/standard_configure.sh` file:
    ```console
    set KRATOS_APPLICATIONS=
    CALL :add_app %KRATOS_APP_DIR%\MPMApplication;
    CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
    ```

## GUI

A GUI (Graphic User Interface) for the MPM application is also available within the pre and post processing software [GiD](https://www.gidhome.com/). Instructions on how to download and install it are available in the [`GiDInterface` repository](https://github.com/KratosMultiphysics/GiDInterface/tree/master/). A basic knowledge of GiD is required.

Any software able to handle `vtk` files can be used for post processing (e.g., [Paraview](https://www.paraview.org/), [VisIt](https://visit-dav.github.io/visit-website/index.html)).

## Examples & Tutorials
* Use-cases and validation examples are available in the MPM section of the [Examples repository](https://kratosmultiphysics.github.io/Examples/).
* Unit tests of the main features can be found in the [tests](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication/tests) folder.
* A step-by-step tutorial using GiD for both pre and post processing is available [here](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Step-by-step_Tutorial_in_GiD/introduction.html).

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
* [Linear isotropic elastic materials](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Constitutive_Laws/constitutive_laws.html#linear-elasticity) - plane strain, plane stress, axis-symmetric, and 3D
* [Hyperelastic Neo-Hookean laws](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Constitutive_Laws/constitutive_laws.html#hyperelastic-neohookean) - finite strain, plane strain, axis-symmetric, and 3D
* Elasto-plastic laws:
    * [Mohr Coulomb](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Constitutive_Laws/constitutive_laws.html#mohr-coulomb) - finite strain, associative and non-associative, plane strain, axis-symmetric, and 3D
    * [Mohr Coulomb with Strain Softening](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Constitutive_Laws/constitutive_laws.html#mohr-coulomb-strain-softening) - finite strain, associative and non-associative, plane strain, axis-symmetric, and 3D
* Critical state laws:
    * [Modified Cam-Clay](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Constitutive_Laws/constitutive_laws.html#modified-cam-clay) - finite strain, plane strain, axis-symmetric, and 3D
    * Johnson Cook Thermal Plastic (just for explicit MPM)
* Displacement-based [Newtonian Fluid Law](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Constitutive_Laws/constitutive_laws.html#newtonian-fluid) - plane strain and 3D

**Boundary conditions**
* Grid-Based Conditions (conforming): applied directly at the background nodes
    * Neumann: [Point load](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Grid-based_Boundary_Conditions/load.html)
    * Neumann: [Line load](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Grid-based_Boundary_Conditions/load.html) (a distributed load applied over a line)
    * Neumann: [Surface load](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Grid-based_Boundary_Conditions/load.html) (a distributed load applied over a face)
    * Dirichlet: [Slip](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Grid-based_Boundary_Conditions/fixed_displacement_boundary_condition.html) and [non-slip](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Grid-based_Boundary_Conditions/slip_boundary_condition.html) conditions for arbitrary inclination
* Material Point-Based Conditions (non-conforming): applied on movable boundary particles
    * Neumann:
        * [moving point load](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Material_Point-based_Boundary_Conditions/point_load.html)
        * interface condition for partitioned coupling with DEM
    * Dirichlet: fixed, slip or contact condition
        * [penalty method](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/Processes/Material_Point-based_Boundary_Conditions/penalty.html)
        * Lagrange multiplier method (*soon in the master branch*)
        * perturbed Lagrangian method (*soon in the master branch*)
        * interface condition for partitioned coupling with FEM, RBS,...

**Time schemes**
* [Implicit](https://kratosmultiphysics.github.io/Kratos/pages/Applications/MPM_Application/MPM_Solver/mpm_implicit_solver.html) - Newmark/Bossak prediction and correction scheme for static, quasi-static, and dynamic problems
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

The `MPMApplication` is **Open Source**. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The **BSD** licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Contact

* **Antonia Larese** - *Group Leader* - [antonia.larese@unipd.it](mailto:antonia.larese@unipd.it)
* **Veronika Singer** - *Developer* - [veronika.singer@tum.de](mailto:veronika.singer@tum.de)
* **Laura Moreno** - *Developer* - [laura.morenomartinez@ua.es](mailto:laura.morenomartinez@ua.es)
* **Andi Makarim Katili** - *Developer* - [andi.katili@tum.de](mailto:andi.katili@tum.de)
* **Nicolò Crescenzio** - *Developer* - [nicolo.crescenzio@math.unipd.it](mailto:nicolo.crescenzio@math.unipd.it)
