---
type: C++ Method
title: ValidateDefaults
resource: sources/kratos_parameters.cpp#L1358-L1403
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_validate_defaults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::ValidateDefaults(const Parameters& rDefaultParameters) const`

# Calls

- [key](../../../../functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_validate_defaults](../../../../functions/tests/test_kratos_parameters/TestParameters/test_validate_defaults.md)