---
type: C++ Method
title: ComputeDistance
resource: modeler/clean_up_problematic_triangles_modeler.cpp#L268-L276
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CleanUpProblematicTrianglesModeler::ComputeDistance( const Node& rNode1, const Node& rNode2 )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CleanUpProblematicGeometries](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries.md)