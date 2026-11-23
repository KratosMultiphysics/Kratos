---
type: C++ Method
title: BuildMasterSlaveConstraints
resource: future/solving_strategies/schemes/implicit_scheme.h#L894-L1011
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
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateConstraintContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/EffectiveEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void BuildMasterSlaveConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData)`

# Calls

- [Start](../../../../../../functions/utilities/timer/Timer/Start.md)
- [BuiltinTimer](../../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [CalculateConstraintContribution](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateConstraintContribution.md)
- [EffectiveEquationId](../../../../../../functions/includes/dof/Dof/EffectiveEquationId.md)
- [Stop](../../../../../../functions/utilities/timer/Timer/Stop.md)