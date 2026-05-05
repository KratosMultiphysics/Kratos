# Structural Mechanics Application

 |             **Application**             |                                                                                    **Description**                                                                                    |                              **Status**                              | **Authors** |
|:---------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------:|:-----------:|
| `StructuralMechanicsApplication` | The *Structural Mechanics Application* contains a series of structural elements, as well as solid elements,  the corresponding strategies, solvers and *Constitutive Laws Application* within *Kratos Multiphysics*. | <img src="https://img.shields.io/badge/Status-%F0%9F%9A%80%20Actively%20developed-Green"  width="300px"> | @KratosMultiphysics/structural-mechanics   |

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/beam_roll_up/data/rollup.gif" alt="Solution" style="width: 300px;"/>
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/use_cases/tensile_test_example/data/animation.gif" alt="Solution" style="width: 300px;"/>
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/beam_shallow_angled_structure/data/shallowAngleBeam.gif" alt="Solution" style="width: 300px;"/>
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/catenoid_formfinding/data/catenoid_normal.gif" alt="Solution" style="width: 300px;"/>
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/four_point_sail_formfinding/data/fourpoint_sail.gif" alt="Solution" style="width: 300px;"/>
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/two_dimensional_circular_truss_arch_snapthrough/data/DispCtrl.gif" alt="Solution" style="width: 300px;"/>
</p>

The application includes tests to check the proper functioning of the application.

## Features:

- **A set of *Neumann* conditions**:
     * *Point loads (loads applied directly on the nodes)*
     * *Point moment (a discret moment applied directly on the nodes)*
     * *Line load (a distributed load applied over a line)*
     * *Surface load (a distributed load applied over a face)*
     * *A simple point contact conditions based on the distance*

- **Solid elements**:
    * *Small displacement elements*
        * Irreducible (pure displacement)
        * Mixed formulation ($BBar$)
        * Mixed formulation ($U-\varepsilon$)
    * *Total Lagrangian elements*
        * Irreducible (pure displacement)
        * Mixed formulation ($U-\Delta V/V$)
        * Mixed formulation ($Q1P0$)
    * *Updated Lagrangian* elements *irreducible (pure displacement)*
    * *Total Lagrangian prismatic solid-shell element (*SPrism*)*

- **Structural elements**:
    * *Zero-dimensional elements* :
        * Nodal concentrated element (both 2D/3D). Includes nodal damping, nodal mass and nodal stiffness
    * *Uni-dimensional elements* :
        * Spring-damper element (3D)
        * Cable element (3D)
        * Truss element (3D)
        * Corrotational beam element (both 2D/3D)
    * *Two-dimensional elements* :
        * Membrane (pre-stressed)
        * Isotropic shell element
        * Thin shell (Quadrilateral and triangular)
        * Thick shell (Quadrilateral and triangular)

- **Constitutive laws**:
    * *Isotropic laws (Plane strain, plane stress and 3D)*
    * *The ones available in [`ConstitutiveLawsApplication`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/ConstitutiveLawsApplication/README.md)*

- **Adjoint Sensitivity Analysis**:
    * *This feature provides the framework to compute sensitivities of structural responses (e.g. displacements, strain energy or stresses) with respect to different types of design variables (e.g. nodal coordinates, material or cross-sectional properties or load intensity) with the adjoint approach*

- **Strategies**:
    * *Formfinding strategies*
    * *Eigensolver strategy*
    * *Harmonic analysis strategies*

- **Schemes**:
    * *Relaxation scheme*
    * *Eigen solver scheme*

- **Convergence criteria**:
    * *For displacement and other *DoF**
    * *For displacement and rotation*

- **Utilities and processes**:
    * *A process to post-process eigenvalues*
    * *A *GiDIO* utility for eigen values*
    * *Process to compute the global mass of the system*
    * *Process to identify the neighbours in a prismatic mesh*
    * *Process to transform a pure shell mesh (local dimension equal to 2), to solid-shell mesh (pure 3D mesh)*

- **+100 Python unittest, including Validation tests, and several cpp tests**

## Examples:

Examples can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/structural_mechanics).