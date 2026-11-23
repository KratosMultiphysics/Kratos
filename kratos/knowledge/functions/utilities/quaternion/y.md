---
type: C++ Function
title: y
resource: utilities/quaternion.h#L128
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/csr_matrix/SpMV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/TransposeSpMV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/LogFittingSlope
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline const T y()const`

# Called by

- [SpMV](../../../functions/containers/csr_matrix/SpMV.md)
- [TransposeSpMV](../../../functions/containers/csr_matrix/TransposeSpMV.md)
- [LogFittingSlope](../../../functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/LogFittingSlope.md)