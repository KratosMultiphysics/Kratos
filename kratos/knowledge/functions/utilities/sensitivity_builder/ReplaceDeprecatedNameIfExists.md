---
type: C++ Function
title: ReplaceDeprecatedNameIfExists
resource: utilities/sensitivity_builder.cpp#L297-L313
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReplaceDeprecatedNameIfExists( Parameters& rSettings, const std::string& rDeprecatedName, const std::string& rNewName)`

# Calls

- [RemoveValue](../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [AddValue](../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [SensitivityBuilder](../../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder.md)