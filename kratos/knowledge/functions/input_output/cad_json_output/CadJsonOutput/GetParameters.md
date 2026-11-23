---
type: C++ Method
title: GetParameters
resource: input_output/cad_json_output.cpp#L20-L41
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CadJsonOutput::GetParameters( ModelPart& rModelPart, Parameters& rCadGeometry, IndexType EchoLevel)`

# Calls

- [AddDouble](../../../../functions/sources/kratos_parameters/Parameters/AddDouble.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [AddEmptyArray](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)