---
type: C++ Method
title: EnforceConsecutiveOrdering
resource: includes/json_io.h#L474-L667
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/json_io/KratosJsonIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EnforceConsecutiveOrdering(rapidjson::Document& d)`

# Calls

- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ReadModelPart](../../../../functions/includes/json_io/KratosJsonIO/ReadModelPart.md)