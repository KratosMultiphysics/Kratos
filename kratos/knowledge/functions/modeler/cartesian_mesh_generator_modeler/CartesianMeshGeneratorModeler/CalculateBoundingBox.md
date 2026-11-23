---
type: C++ Method
title: CalculateBoundingBox
resource: modeler/cartesian_mesh_generator_modeler.h#L390-L421
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/point/Point/Dimension
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateBoundingBox(ModelPart& rThisModelPart, Point& rMinPoint, Point& rMaxPoint)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Dimension](../../../../functions/geometries/point/Point/Dimension.md)