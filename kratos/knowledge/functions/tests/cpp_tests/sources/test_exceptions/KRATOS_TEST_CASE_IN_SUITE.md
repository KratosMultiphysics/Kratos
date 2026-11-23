---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/sources/test_exceptions.cpp#L28-L39
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/exception/Exception/what
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/exception/Exception/where
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/CleanFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/CleanFunctionName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/GetLineNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ExceptionDefaultConstruction, KratosCoreFastSuite)`

# Calls

- [what](../../../../../functions/sources/exception/Exception/what.md)
- [where](../../../../../functions/sources/exception/Exception/where.md)
- [CleanFileName](../../../../../functions/sources/code_location/CodeLocation/CleanFileName.md)
- [CleanFunctionName](../../../../../functions/sources/code_location/CodeLocation/CleanFunctionName.md)
- [GetLineNumber](../../../../../functions/sources/code_location/CodeLocation/GetLineNumber.md)