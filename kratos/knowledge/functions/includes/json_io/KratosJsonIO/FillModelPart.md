---
type: C++ Method
title: FillModelPart
resource: includes/json_io.h#L255-L468
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
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
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

`void FillModelPart(rapidjson::Document& d, ModelPart & rThisModelPart)`

# Calls

- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ReadModelPart](../../../../functions/includes/json_io/KratosJsonIO/ReadModelPart.md)