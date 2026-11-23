---
type: C++ Method
title: AddStringArray
resource: sources/kratos_parameters.cpp#L961-L969
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_string_array_valid
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddStringArray( const std::string& rEntry, const std::vector<std::string>& rValue )`

# Calls

- [SetStringArray](../../../../functions/sources/kratos_parameters/Parameters/SetStringArray.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [test_add_string_array_valid](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_string_array_valid.md)