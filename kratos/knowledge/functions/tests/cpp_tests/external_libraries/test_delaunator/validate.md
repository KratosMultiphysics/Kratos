---
type: C++ Function
title: validate
resource: tests/cpp_tests/external_libraries/test_delaunator.cpp#L100-L107
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void validate( const std::vector<double>& rCoordinates, const double HullArea = 0.0, const std::string& Debug = "" )`

# Calls

- [validatewithtolerance](../../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/KRATOS_TEST_CASE_IN_SUITE.md)