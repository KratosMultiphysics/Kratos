---
type: C++ Method
title: AddDouble
resource: sources/kratos_parameters.cpp#L909-L917
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_methods
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddDouble( const std::string& rEntry, const double Value )`

# Calls

- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [GetParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetParameters.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_add_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_methods.md)