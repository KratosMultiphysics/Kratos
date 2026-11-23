---
type: C++ Method
title: GetLineNumber
resource: sources/code_location.cpp#L37-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/input_output/test_logger/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/profiler/Profiler/Write
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int CodeLocation::GetLineNumber() const`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/input_output/test_logger/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE.md)
- [Write](../../../../functions/utilities/profiler/Profiler/Write.md)