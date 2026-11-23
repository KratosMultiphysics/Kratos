---
type: C++ Function
title: TrimRight
resource: utilities/string_utilities.cpp#L227-L240
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/string_utilities/TrimChar
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/Trim
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string TrimRight( const std::string& rInputString, const bool RemoveNullChar)`

# Calls

- [TrimChar](../../../functions/utilities/string_utilities/TrimChar.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_string_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [Trim](../../../functions/utilities/string_utilities/Trim.md)