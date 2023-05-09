# Kratos ROM Application

The ROM application in Kratos Multiphysics is a powerful tool for reducing the computational resources required to simulate complex engineering problems. It includes advanced mathematical techniques such as Proper Orthogonal Decomposition (POD), Empirical Cubature Method, Randomized SVD, and MPI SVD to create reduced-order models that accurately simulate the behavior of a system with a significantly lower number of degrees of freedom than the original system.

<img src="https://github.com/KratosMultiphysics/Kratos/assets/61457043/42ad4f92-6fdd-4509-bbd6-ee0f61a1592d" alt="ROM_logo" width="200"/>

## Rom Manager

The Kratos ROM application comes equipped with a Rom Manager, which provides users with the ability to reduce their models using a variety of techniques. These include reducing the model using the Proper Orthogonal Decomposition (POD) technique with Galerkin, Least-Squares Petrov Galerkin, and Petrov-Galerkin projections, as well as hyper-reducing the model using the Empirical Cubature Method (ECM). The Rom Manager is available in the Kratos Multiphysics Examples repository [link](https://github.com/KratosMultiphysics/Examples/tree/master/rom_application/RomManager).

## Compatibility

The ROM application is compatible with and has been tested on the following applications in Kratos Multiphysics:

- FluidDynamicsApplication
- StructuralMechanicsApplication
- ConvectionDiffusionApplication
- CompressiblePotentialFlowApplication
- CoSimulationApplication

## Applications

The ROM application has various applications in several fields, including:

- Fluid dynamics
- Structural mechanics
- Electromagnetics

## Planned Expansions

The ROM application is designed to be expanded to a parallel High-Performance Computing (HPC) workflow, allowing for even more efficient computation of reduced-order models on large-scale systems. Additionally, the ROM application is planned to be able to use local basis, further increasing the accuracy and efficiency of the reduced-order models.

## Availability

The ROM application is available on [GitHub](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/RomApplication) for researchers and engineers to use and contribute to.
