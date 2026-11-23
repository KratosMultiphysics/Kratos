---
type: C++ Method
title: SetupGeometryModel
resource: modeler/cad_io_modeler.cpp#L24-L68
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
---

# Signature

`void CadIoModeler::SetupGeometryModel()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)