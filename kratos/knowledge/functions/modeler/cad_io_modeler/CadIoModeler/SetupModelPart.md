---
type: C++ Method
title: SetupModelPart
resource: modeler/cad_io_modeler.cpp#L70-L87
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/GetCadJsonOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CadIoModeler::SetupModelPart()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [GetCadJsonOutput](../../../../functions/input_output/cad_json_output/GetCadJsonOutput.md)