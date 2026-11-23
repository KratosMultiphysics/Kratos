---
type: C++ Method
title: NodesOfBoundaries
resource: processes/tetrahedral_mesh_orientation_check.cpp#L490-L514
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/BoundariesEntitiesNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NumberOfNodesInEachBoundary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/LinearElement
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TetrahedralMeshOrientationCheck::NodesOfBoundaries(const GeometryType& rGeometry, DenseMatrix<int>& rNodesIds)`

# Calls

- [BoundariesEntitiesNumber](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/BoundariesEntitiesNumber.md)
- [NumberOfNodesInEachBoundary](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NumberOfNodesInEachBoundary.md)
- [LinearElement](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/LinearElement.md)

# Called by

- [Execute](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute.md)