---
type: C++ Method
title: GenerateCutGeometryNodes
resource: modeler/entity_generation/generate_stl_intersection_with_cells.cpp#L89-L102
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`PointerVector<Node> GenerateStlIntersectionWithCells::GenerateCutGeometryNodes( ModelPart& rTheCutModelPart, const std::vector<typename WorkGeometryType::Pointer>& rCutGeometries)`

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate.md)