---
type: C++ Method
title: BuildRHSAndSolve
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L660-L693
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BuildRHSAndSolve( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)
- [BuiltinTimer](../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)