---
type: C++ Method
title: GenerateIntersection
resource: modeler/entity_generation/generate_stl_intersection_with_cells.cpp#L111-L135
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/skin_intersection/ConvertToWorkGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/skin_intersection/SkinIntersection/ComputeCutInsideCell
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateStlIntersectionWithCells::GenerateIntersection(const ModelPart& rReferenceModelPart, const std::string& rInputEntities, int Color)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ConvertToWorkGeometry](../../../../../functions/modeler/internals/skin_intersection/ConvertToWorkGeometry.md)
- [CalculateOuterMinMaxNodePositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [ComputeCutInsideCell](../../../../../functions/modeler/internals/skin_intersection/SkinIntersection/ComputeCutInsideCell.md)