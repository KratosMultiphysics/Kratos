---
type: C++ Method
title: PerformRungeKuttaLastSubStep
resource: solving_strategies/strategies/explicit_solving_strategy_runge_kutta.h#L417-L450
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/SubstepCount
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetIntegrationTheta
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void PerformRungeKuttaLastSubStep(LocalSystemMatrixType& rLastStepResidualVector)`

# Calls

- [SubstepCount](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/SubstepCount.md)
- [GetIntegrationTheta](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetIntegrationTheta.md)

# Called by

- [SolveWithLumpedMassMatrix](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix.md)