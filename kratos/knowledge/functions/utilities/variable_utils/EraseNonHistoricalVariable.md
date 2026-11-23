---
type: C++ Function
title: EraseNonHistoricalVariable
resource: utilities/variable_utils.h#L838-L850
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_variable_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EraseNonHistoricalVariable( const TVarType& rVariable, TContainerType& rContainer )`

# Calls

- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_variable_utils/KRATOS_TEST_CASE_IN_SUITE.md)