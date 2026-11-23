---
type: C++ Method
title: NumberOfNodesInEachBoundary
resource: processes/tetrahedral_mesh_orientation_check.cpp#L472-L485
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NodesOfBoundaries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`TetrahedralMeshOrientationCheck::SizeType TetrahedralMeshOrientationCheck::NumberOfNodesInEachBoundary(const GeometryType& rGeometry)`

# Called by

- [Execute](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute.md)
- [NodesOfBoundaries](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NodesOfBoundaries.md)