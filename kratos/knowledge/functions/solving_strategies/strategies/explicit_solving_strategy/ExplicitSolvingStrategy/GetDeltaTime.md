---
type: C++ Method
title: GetDeltaTime
resource: solving_strategies/strategies/explicit_solving_strategy.h#L410-L413
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/SolveWithLumpedMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/PerformSubstep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual inline double GetDeltaTime()`

# Called by

- [SolveWithLumpedMassMatrix](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/SolveWithLumpedMassMatrix.md)
- [PerformSubstep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/PerformSubstep.md)
- [SolveWithLumpedMassMatrix](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix.md)
- [PerformRungeKuttaIntermediateSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep.md)