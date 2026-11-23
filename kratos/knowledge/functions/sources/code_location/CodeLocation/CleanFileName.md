---
type: C++ Method
title: CleanFileName
resource: sources/code_location.cpp#L41-L57
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/code_location/CodeLocation/ReplaceAll
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_code_location/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string CodeLocation::CleanFileName() const`

# Calls

- [ReplaceAll](../../../../functions/sources/code_location/CodeLocation/ReplaceAll.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_code_location/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE.md)