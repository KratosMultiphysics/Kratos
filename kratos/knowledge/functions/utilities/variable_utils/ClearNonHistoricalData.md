---
type: C++ Function
title: ClearNonHistoricalData
resource: utilities/variable_utils.h#L857-L866
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
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_clear_nonhistorical_values
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ClearNonHistoricalData(TContainerType& rContainer)`

# Calls

- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_variable_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [setUp](../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/setUp.md)
- [setUp](../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/setUp.md)
- [test_clear_nonhistorical_values](../../../functions/tests/test_variable_utils/TestVariableUtils/test_clear_nonhistorical_values.md)