---
type: C++ Method
title: BuildAndSolve
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h#L508-L541
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BuildAndSolve( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb ) override`

# Calls

- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)
- [BuiltinTimer](../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [Start](../../../../../functions/utilities/timer/Timer/Start.md)