---
type: C++ Module
title: structured_mesh_generator_process
resource: processes/structured_mesh_generator_process.cpp#L1-L421
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/vector
    resolved_by: tree-sitter
    confidence: exact
  - target: external/numeric
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-structured-mesh-generator-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-point-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-checks-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-skin-detection-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-interface-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [StructuredMeshGeneratorProcess](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/StructuredMeshGeneratorProcess.md)
- [Execute](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Execute.md)
- [ValidateTheDefaultParameters](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters.md)
- [GetDefaultParameters](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetDefaultParameters.md)
- [Info](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Info.md)
- [PrintInfo](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/PrintInfo.md)
- [PrintData](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/PrintData.md)
- [Generate2DMesh](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh.md)
- [Generate3DMesh](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh.md)
- [GenerateNodes2D](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateNodes2D.md)
- [GenerateNodes3D](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateNodes3D.md)
- [GenerateTriangularElements](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements.md)
- [GenerateTetrahedraElements](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements.md)
- [CreateCellTetrahedra](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CreateCellTetrahedra.md)
- [GetNodeId](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetNodeId.md)
- [GetLocalCoordinatesRange](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetLocalCoordinatesRange.md)
- [Check](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Check.md)
- [CheckDomainGeometry](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CheckDomainGeometry.md)
- [CheckDomainGeometryConnectivityForQuadrilateral2D4](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CheckDomainGeometryConnectivityForQuadrilateral2D4.md)
- [CheckDomainGeometryConnectivityForHexahedra3D8](../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CheckDomainGeometryConnectivityForHexahedra3D8.md)

# Imports

- `<vector>`
- `<numeric>`
- `"processes/structured_mesh_generator_process.h"`
- `"geometries/geometry.h"`
- `"geometries/point.h"`
- `"geometries/triangle_2d_3.h"`
- `"geometries/tetrahedra_3d_4.h"`
- `"includes/checks.h"`
- `"processes/skin_detection_process.h"`
- `"includes/ublas_interface.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)