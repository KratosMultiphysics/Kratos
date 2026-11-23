---
type: C++ Method
title: BuildDampingMatrix
resource: future/solving_strategies/schemes/implicit_scheme.h#L833-L892
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
  - target: functions/future/solving_strategies/schemes/implicit_scheme/CalculateDampingMatrixContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void BuildDampingMatrix(MatrixType& rDampingMatrix)`

# Calls

- [Start](../../../../../../functions/utilities/timer/Timer/Start.md)
- [BuiltinTimer](../../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [CalculateDampingMatrixContribution](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/CalculateDampingMatrixContribution.md)
- [Stop](../../../../../../functions/utilities/timer/Timer/Stop.md)