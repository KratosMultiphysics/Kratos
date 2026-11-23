---
type: C++ Method
title: ReadControlPointVector
resource: input_output/cad_json_input.h#L798-L814
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadAndCreateNode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadControlPointVector( PointerVector<Node>& rControlPoints, const Parameters rParameters, ModelPart& rModelPart, SizeType EchoLevel = 0)`

# Calls

- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [ReadAndCreateNode](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadAndCreateNode.md)