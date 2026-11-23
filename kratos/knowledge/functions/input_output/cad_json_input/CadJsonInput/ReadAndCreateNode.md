---
type: C++ Method
title: ReadAndCreateNode
resource: input_output/cad_json_input.h#L843-L857
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadControlPointVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static Node::Pointer ReadAndCreateNode( const Parameters rParameters, ModelPart& rModelPart, SizeType EchoLevel = 0)`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)

# Called by

- [ReadControlPointVector](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadControlPointVector.md)