---
type: C++ Method
title: CopyValuesFromExistingParameters
resource: sources/kratos_parameters.cpp#L1198-L1214
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_kratos_parameters/TestParameters/test_copy_values_from_existing_parameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::CopyValuesFromExistingParameters( const Parameters OriginParameters, const std::vector<std::string>& rListParametersToCopy )`

# Calls

- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [test_copy_values_from_existing_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_copy_values_from_existing_parameters.md)