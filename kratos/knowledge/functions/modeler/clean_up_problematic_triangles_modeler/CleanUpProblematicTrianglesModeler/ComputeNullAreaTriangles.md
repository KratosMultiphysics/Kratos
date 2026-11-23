---
type: C++ Method
title: ComputeNullAreaTriangles
resource: modeler/clean_up_problematic_triangles_modeler.cpp#L298-L318
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeSquaredArea
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t CleanUpProblematicTrianglesModeler::ComputeNullAreaTriangles( ModelPart& rThisModelPart, const double RefArea, const double AreaTolerance )`

# Calls

- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [ComputeSquaredArea](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeSquaredArea.md)