---
type: C++ Method
title: GetMatrixRowBegin
resource: solving_strategies/strategies/butcher_tableau.h#L122-L126
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRow
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`RowType::const_iterator GetMatrixRowBegin(const unsigned int SubStepIndex) const`

# Called by

- [GetMatrixRow](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRow.md)
- [PerformRungeKuttaIntermediateSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep.md)