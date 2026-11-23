---
type: C++ Method
title: Execute
resource: processes/tetrahedral_mesh_orientation_check.cpp#L142-L358
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SupportedElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Orient
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SupportedCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/BoundariesEntitiesNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NumberOfNodesInEachBoundary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NodesOfBoundaries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TetrahedralMeshOrientationCheck::Execute()`

# Calls

- [SolutionStepsDataHas](../../../../functions/includes/node/SolutionStepsDataHas.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SupportedElement](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SupportedElement.md)
- [Orient](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Orient.md)
- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [SupportedCondition](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SupportedCondition.md)
- [BoundariesEntitiesNumber](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/BoundariesEntitiesNumber.md)
- [NumberOfNodesInEachBoundary](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NumberOfNodesInEachBoundary.md)
- [NodesOfBoundaries](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/NodesOfBoundaries.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)