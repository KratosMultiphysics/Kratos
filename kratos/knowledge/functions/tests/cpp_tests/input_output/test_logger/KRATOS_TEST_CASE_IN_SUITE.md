---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/input_output/test_logger.cpp#L32-L58
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/logger_message/GetLocation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/GetLineNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_output/SetOption
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(LoggerMessageStream, KratosCoreFastSuite)`

# Calls

- [GetLocation](../../../../../functions/input_output/logger_message/GetLocation.md)
- [GetLineNumber](../../../../../functions/sources/code_location/CodeLocation/GetLineNumber.md)
- [SetOption](../../../../../functions/input_output/logger_output/SetOption.md)