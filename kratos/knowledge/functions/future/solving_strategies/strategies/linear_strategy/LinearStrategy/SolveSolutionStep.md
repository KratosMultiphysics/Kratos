---
type: C++ Method
title: SolveSolutionStep
resource: future/solving_strategies/strategies/linear_strategy.h#L183-L263
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/pGetScheme
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/pGetLinearSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/GetComputeReactions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SolveSolutionStep() override`

# Calls

- [pGetScheme](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/pGetScheme.md)
- [pGetVector](../../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [BuildLinearSystemConstraints](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildLinearSystemConstraints.md)
- [pGetLinearSolver](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/pGetLinearSolver.md)
- [GetComputeReactions](../../../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/GetComputeReactions.md)