---
type: C++ Method
title: RemoveValues
resource: sources/kratos_parameters.cpp#L490-L499
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_kratos_parameters/TestParameters/test_remove_values
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::RemoveValues(const std::vector<std::string>& rEntries)`

# Calls

- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)

# Called by

- [test_remove_values](../../../../functions/tests/test_kratos_parameters/TestParameters/test_remove_values.md)