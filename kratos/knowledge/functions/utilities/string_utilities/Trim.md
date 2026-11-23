---
type: C++ Function
title: Trim
resource: utilities/string_utilities.cpp#L182-L187
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/TrimLeft
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/TrimRight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/StringToVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string Trim( const std::string& rInputString, const bool RemoveNullChar)`

# Calls

- [TrimLeft](../../../functions/utilities/string_utilities/TrimLeft.md)
- [TrimRight](../../../functions/utilities/string_utilities/TrimRight.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [StringToVector](../../../functions/utilities/string_utilities/StringToVector.md)