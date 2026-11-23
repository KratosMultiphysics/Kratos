---
type: C++ Method
title: BuildAndSolveLinearizedOnPreviousIteration
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L552-L650
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FreeDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/VariableUtils/UpdateCurrentPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BuildAndSolveLinearizedOnPreviousIteration( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb, const bool MoveMesh ) override`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [FreeDof](../../../../../functions/includes/dof/Dof/FreeDof.md)
- [UpdateCurrentPosition](../../../../../functions/utilities/variable_utils/VariableUtils/UpdateCurrentPosition.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)
- [FixDof](../../../../../functions/includes/dof/Dof/FixDof.md)
- [BuiltinTimer](../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)