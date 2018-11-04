# Meshing Application

Meshing Application provides several tools to create, manipulate and interact with meshes. It containts several interfaces to both Kratos
thrid party libraries (Triangle, TetGen, MMG)

The application offers the functionalities listed below. If there is an Object without methds it means it can be called using the __Execute()__ function.

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
  
## Interface

### Custom IO
* __PFEMGidIO__: A specialized instance of GiDIO for the PFEM Application. It redefines several IO methods:
  * _WriteMesh_
  * _WriteNodeMesh_
  * _InitializeMesh_
  * _FinalizeMesh_
  * _InitializeResults_
  * _FinalizeResults_
  * _WriteNodalResults_
  * _PrintOnGaussPoints_
  * _Flush_
  * _CloseResultFile_
  
### Utilities
* __MeshTransfer2D__
* __MeshTransfer3D__:
  * _DirectModelPartInterpolation_
  * _DirectScalarVarInterpolation_
  * _DirectVectorialVarInterpolation_

* __BinBasedMeshTransfer2D__
* __BinBasedMeshTransfer3D__: Alternative implementation of the __MeshTransfer__ utility based on bins. Inherits the procedures from __MeshTransfer__ and also adds:
  * _MappingFromMovingMesh_ScalarVar_
  * _MappingFromMovingMesh_VectorialVar_
  * _MappingFromMovingMesh_VariableMeshes_ScalarVar_
  * _MappingFromMovingMesh_VariableMeshes_VectorialVar_

* __LocalRefineTriangleMesh__: Refines a Trianglular Mesh.
* __LocalRefinePrismMesh__: Refines a Prism Mesh.
* __LocalRefineSPrismMesh__: Refines a SPrism Mesh.
* __LocalRefineTetrahedraMesh__: Refines a Tetrahedra Mesh.

* __Cutting_Isosurface_Application__:
  * _GenerateScalarVarCut_
  * _GenerateVectorialComponentVarCut_
  * _GenerateVectorialVarCut_
  * _AddModelPartElements_
  * _AddSkinConditions_
  * _UpdateCutData_
  * _DeleteCutData_

### Meshers
* __TriGenPFEMModeler__:
  * _ReGenerateMesh_
  
* __TriGenGLASSModeler__:
  * _ReGenerateMeshGlass_
  
* __TriGenPFEMModelerVMS__:
  * _ReGenerateMesh_
  
* __TriGenPFEMSegment__:
  * _ReGenerateMesh_

### Processes
* __InternalVariablesInterpolationProcess__: Inerpolates Nodal v

#### Metrics
* __MetricFastInit2D__:
* __MetricFastInit3D__:

#### LevelSet
* __ComputeLevelSetSolMetricProcess2D__:
* __ComputeLevelSetSolMetricProcess3D__:
        
#### Hessian
* __ComputeHessianSolMetricProcess2D__: For double values.
* __ComputeHessianSolMetricProcess3D__: For double values.
* __ComputeHessianSolMetricProcessComp2D__: For components.
* __ComputeHessianSolMetricProcessComp3D__: For components.
        
#### Error
* __ComputeErrorSolMetricProcess2D__:
* __ComputeErrorSolMetricProcess3D__:

## External Libraries

Meshing application can make use of several third party libs as an alternative (or sometimes unique) way to implementd the
interface shown. You can find information about these libs in their respective pages which are listed below:

### TetGen
[TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) is a program to generate tetrahedral meshes of any 3D polyhedral domains. 
Please note that __Tetgen license is not compatible__ with Kratos, and hence it is not included as part of Kratos. You must indicate in compile time where it can find a tetgen already in your system.

Tetgen enables to use the following __utilities__:

* __TetgenVolumeMesher__:
  * _AddHole_
  * _GenerateMesh_
    
* __TetrahedraReconnectUtility__:
  * _EvaluateQuality_
  * _TestRemovingElements_
  * _OptimizeQuality_
  * _FinalizeOptimization_
  * _updateNodesPositions_
  * _setMaxNumThreads_
  * _setBlockSize_
  * _isaValidMesh_
  
Tetgen also enable to use the following __meshers__:

* __TetGenPfemModeler__:
  * _ReGenerateMesh_

* __TetGenPfemRefineFace__:
  * _ReGenerateMesh_

* __TetGenPfemContact__:
  * _ReGenerateMesh_
    
* __TetGenCDT__:
  * _GenerateCDT_

* __TetGenPfemModelerVms__:
  * _ReGenerateMesh_
  
### MMG
[MMG](https://www.mmgtools.org/) is an open source software for simplicial remeshing. It provides 3 applications and 4 libraries. In Kratos it provides the following additioanl __procedures__:

* __MmgProcess2D__:
* __MmgProcess3D__:
        
