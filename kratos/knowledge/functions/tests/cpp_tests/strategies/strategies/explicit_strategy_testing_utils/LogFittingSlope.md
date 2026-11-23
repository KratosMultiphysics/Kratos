---
type: C++ Function
title: LogFittingSlope
resource: tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils.h#L170-L190
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/quaternion/y
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline double LogFittingSlope( const std::array<double, TDataSize>& rX, const std::array<double, TDataSize>& rY)`

# Calls

- [y](../../../../../../functions/utilities/quaternion/y.md)
- [prod](../../../../../../functions/includes/amatrix_interface/prod.md)