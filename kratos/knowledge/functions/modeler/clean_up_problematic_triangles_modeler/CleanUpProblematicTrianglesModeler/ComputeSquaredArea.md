---
type: C++ Method
title: ComputeSquaredArea
resource: modeler/clean_up_problematic_triangles_modeler.cpp#L281-L292
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/Norm3
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeNullAreaTriangles
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CleanUpProblematicTrianglesModeler::ComputeSquaredArea(const GeometryType& rGeometry)`

# Calls

- [Norm3](../../../../functions/utilities/math_utils/Norm3.md)

# Called by

- [CleanUpProblematicGeometries](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries.md)
- [ComputeNullAreaTriangles](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeNullAreaTriangles.md)