# Meshing Application

 |            **Application**             |                                                                                    **Description**                                                                                    |                              **Status**                              | **Authors** |
|:---------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------:|:-----------:|
| `MeshingApplication` | The *Meshing Application* is a core component of the Kratos Multiphysics framework that provides tools and processes for mesh generation, adaptation, refinement, and interact with meshes within *Kratos Multiphysics*. | <img src="https://img.shields.io/badge/Status-%F0%9F%9A%80%20Actively%20developed-Green"  width="300px"> | [*Pooyan Dadvand*](mailto:pooyan@altair.com)  <br />  [*Vicente Mataix Ferrándiz*](mailto:vmataix@altair.com) |

## Features

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/blob/master/mmg_remeshing_examples/use_cases/channel_sphere2D/data/result.gif?raw=true" alt="Solution" style="width: 600px;"/>
</p>

It offers various algorithms for computing metrics, remeshing based on different criteria, and interpolating variables between meshes.

It contains several interfaces to both Kratos third party libraries (*Triangle*, *TetGen*, *MMG*)

The application offers the functionalities listed below. If there is an Object without methods it means it can be called using the `Execute()` function.

- [Interface](#interface)
  * [Custom IO](#custom-io)
  * [Utilities](#utilities)
  * [Meshers](#meshers)
  * [Processes](#processes)
    + [Metrics](#metrics)
    + [LevelSet](#levelset)
    + [Hessian](#hessian)
    + [Error](#error)
- [External Libraries](#external-libraries)
  * [TetGen](#tetgen)
  * [MMG](#mmg)

### Custom IO
* `PFEMGidIO`: A specialized instance of GiDIO for the PFEM Application. 
  
### Utilities
* `MeshTransfer2D/3D`
* `BinBasedMeshTransfer2D`: Alternative implementation of the `MeshTransfer` utility based on bins. 
* `LocalRefineTriangleMesh`: Refines a Triangular Mesh.
* `LocalRefinePrismMesh`: Refines a Prism Mesh.
* `LocalRefineSPrismMesh`: Refines a SPrism Mesh.
* `LocalRefineTetrahedraMesh`: Refines a Tetrahedra Mesh.
* `Cutting_Isosurface_Application`

### Meshers
* `TriGenPFEMModeler`
* `TriGenGLASSModeler`
* `TriGenPFEMModelerVMS`
* `TriGenPFEMSegment`

### Processes

#### Variable Interpolation

After remeshing, variables need to be transferred from the old mesh to the new one. The application provides two main interpolation processes:

1. `NodalValuesInterpolationProcess`: Interpolates nodal solution step variables.
2. `InternalVariablesInterpolationProcess`: Interpolates internal variables stored at integration points. Options are:
  - `CLOSEST_POINT_TRANSFER`
  - `LEAST_SQUARE_TRANSFER`
  - `SHAPE_FUNCTION_TRANSFER`

#### Metrics

Metrics are represented as tensors that define the desired mesh size and direction at each point. For 2D problems, a 3-component symmetric tensor is used, while for 3D problems, a 6-component symmetric tensor is used.

* `MetricFastInit2D/3D`: Initializes the metric tensors.
1. LevelSet
  * `ComputeLevelSetSolMetricProcess2D/3D`: Computes the level-set related metric.
2. Hessian: Computes the *Hessian* based metric.
  * `ComputeHessianSolMetricProcess2D/3D`: For double values.
  * `ComputeHessianSolMetricProcessComp2D/3D`: For components.
3. Error
  * `ComputeErrorSolMetricProcess2D/3D`: Computes the metric associated to the error computed using **SPR** (*Super Patch Recovery*).

### External Libraries

Meshing application can make use of several third party libs as an alternative (or sometimes unique) way to implemented the interface shown. You can find information about these libs in their respective pages which are listed below:

#### TetGen

[TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) is a program to generate tetrahedral meshes of any 3D polyhedral domains. 
Please note that `Tetgen license is not compatible` with Kratos, and hence it is not included as part of Kratos. You must indicate in compile time where it can find a tetgen already in your system.

*Tetgen* enables to use the following `utilities`:

* `TetgenVolumeMesher`
* `TetrahedraReconnectUtility`
  
*Tetgen* also enable to use the following `meshers`:

* `TetGenPfemModeler`
* `TetGenPfemRefineFace`
* `TetGenPfemContact`
* `TetGenCDT`
* `TetGenPfemModelerVms`
  
#### MMG

[MMG](https://www.mmgtools.org/) is an open source software for simplicial remeshing. It provides 3 applications and 4 libraries. In Kratos it provides the following additional `procedures`:

* `MmgProcess`: This class is a remesher which uses the MMG library. The class uses a class for the 2D and 3D cases (solid and surfaces). The remesher keeps the previous submodelparts and interpolates the nodal values between the old and new mesh.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/MMG_files/mmg_workflow.png?raw=true" alt="MMG architecture" style="width: 600px;"/>
</p>

##### Architecture

The following diagram illustrates the high-level architecture of the *MMG* remeshing:

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/MMG_files/mmg_architecture.png?raw=true" alt="MMG architecture" style="width: 600px;"/>
</p>

##### Remeshing Strategies

The application supports several remeshing strategies:

1. *Hessian-based*: Adapts the mesh based on the Hessian of a scalar field, which represents the curvature of the solution.
2. *Level-set*: Adapts the mesh based on a distance function gradient, often used for interface problems.
3. *Error-based*: Adapts the mesh based on error estimators or indicators.
4. *Optimization*: Only optimizes the mesh quality without significant adaptation.

## Examples:

### MMG examples

Examples can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/mmg_remeshing_examples).

### ParMMG examples

Examples can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/parmmg_remeshing_examples).


