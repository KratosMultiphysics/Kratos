---
type: C++ Function
title: SetReshapeMatrixFlag
resource: solving_strategies/builder_and_solvers/builder_and_solver.h#L229-L232
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/AdaptiveResidualBasedNewtonRaphsonStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/SetReformDofSetAtEachStepFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/ResidualBasedLinearStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/SetReformDofSetAtEachStepFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/ResidualBasedNewtonRaphsonStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SetReformDofSetAtEachStepFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetReshapeMatrixFlag(bool ReshapeMatrixFlag)`

# Called by

- [AdaptiveResidualBasedNewtonRaphsonStrategy](../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/AdaptiveResidualBasedNewtonRaphsonStrategy.md)
- [SetReformDofSetAtEachStepFlag](../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/SetReformDofSetAtEachStepFlag.md)
- [ResidualBasedLinearStrategy](../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/ResidualBasedLinearStrategy.md)
- [SetReformDofSetAtEachStepFlag](../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/SetReformDofSetAtEachStepFlag.md)
- [ResidualBasedNewtonRaphsonStrategy](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/ResidualBasedNewtonRaphsonStrategy.md)
- [SetReformDofSetAtEachStepFlag](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SetReformDofSetAtEachStepFlag.md)