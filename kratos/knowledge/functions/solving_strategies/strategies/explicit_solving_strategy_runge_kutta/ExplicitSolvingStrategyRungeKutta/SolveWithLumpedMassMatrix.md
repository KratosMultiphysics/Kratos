---
type: C++ Method
title: SolveWithLumpedMassMatrix
resource: solving_strategies/strategies/explicit_solving_strategy_runge_kutta.h#L257-L315
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/SubstepCount
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaLastSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SolveWithLumpedMassMatrix() override`

# Calls

- [GetDeltaTime](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime.md)
- [SubstepCount](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/SubstepCount.md)
- [PerformRungeKuttaIntermediateSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep.md)
- [PerformRungeKuttaLastSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaLastSubStep.md)