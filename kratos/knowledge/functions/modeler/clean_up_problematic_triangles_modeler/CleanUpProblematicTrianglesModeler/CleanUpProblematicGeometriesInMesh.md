---
type: C++ Method
title: CleanUpProblematicGeometriesInMesh
resource: modeler/clean_up_problematic_triangles_modeler.cpp#L49-L67
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CleanUpProblematicTrianglesModeler::CleanUpProblematicGeometriesInMesh( ModelPart& rThisModelPart, const std::string& rEntityType, const IndexType FirstNodeId, const IndexType FirstElementId, const IndexType FirstConditionId, const double AreaTolerance )`

# Called by

- [SetupModelPart](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/SetupModelPart.md)