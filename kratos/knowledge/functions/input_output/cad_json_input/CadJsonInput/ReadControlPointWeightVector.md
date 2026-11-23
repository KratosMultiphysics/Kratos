---
type: C++ Method
title: ReadControlPointWeightVector
resource: input_output/cad_json_input.h#L775-L795
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/ReadNurbsCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/ReadNurbsSurface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static Vector ReadControlPointWeightVector( const Parameters rParameters, SizeType EchoLevel = 0)`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)

# Called by

- [ReadNurbsCurve](../../../../functions/input_output/cad_json_input/ReadNurbsCurve.md)
- [ReadNurbsSurface](../../../../functions/input_output/cad_json_input/ReadNurbsSurface.md)