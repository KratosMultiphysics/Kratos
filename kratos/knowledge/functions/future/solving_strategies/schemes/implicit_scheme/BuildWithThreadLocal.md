---
type: C++ Function
title: BuildWithThreadLocal
resource: future/solving_strategies/schemes/implicit_scheme.h#L484-L566
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/tests/cpp_tests/test_schemes/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void BuildWithThreadLocal( MatrixType& rLHS, VectorType& rRHS)`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [BuiltinTimer](../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/future/tests/cpp_tests/test_schemes/KRATOS_TEST_CASE_IN_SUITE.md)