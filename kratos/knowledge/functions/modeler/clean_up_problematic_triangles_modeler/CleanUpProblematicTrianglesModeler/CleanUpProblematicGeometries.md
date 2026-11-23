---
type: C++ Method
title: CleanUpProblematicGeometries
resource: modeler/clean_up_problematic_triangles_modeler.cpp#L73-L238
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeSquaredArea
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CleanUpProblematicTrianglesModeler::CleanUpProblematicGeometries( ModelPart& rThisModelPart, const IndexType FirstNodeId, const IndexType FirstEntityId, const double AreaTolerance )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [ComputeSquaredArea](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeSquaredArea.md)
- [ComputeDistance](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeDistance.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)