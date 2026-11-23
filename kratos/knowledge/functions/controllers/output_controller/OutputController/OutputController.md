---
type: C++ Method
title: OutputController
resource: controllers/output_controller.cpp#L34-L64
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`OutputController::OutputController( const Model& rModel, Parameters Settings) : mpModel(&rModel)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)