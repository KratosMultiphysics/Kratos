---
type: C++ Method
title: BuildMassMatrix
resource: future/solving_strategies/schemes/implicit_scheme.h#L772-L831
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
  - target: functions/future/solving_strategies/schemes/implicit_scheme/CalculateMassMatrixContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void BuildMassMatrix(MatrixType& rMassMatrix)`

# Calls

- [Start](../../../../../../functions/utilities/timer/Timer/Start.md)
- [BuiltinTimer](../../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [CalculateMassMatrixContribution](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateMassMatrixContribution.md)
- [Stop](../../../../../../functions/utilities/timer/Timer/Stop.md)