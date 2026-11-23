---
type: C++ Function
title: CheckAndAssignNodalValue
resource: includes/json_io.h#L670-L696
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CheckAndAssignNodalValue( ModelPart& r_model_part, std::string variable_name, const rapidjson::Value& json_value, unsigned int stepindex)`

# Calls

- [GetInt](../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [Fix](../../../functions/includes/node/Fix.md)