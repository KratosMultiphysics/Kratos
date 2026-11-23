---
type: C++ Method
title: PerformRungeKuttaIntermediateSubStep
resource: solving_strategies/strategies/explicit_solving_strategy_runge_kutta.h#L349-L410
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetIntegrationTheta
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowEnd
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void PerformRungeKuttaIntermediateSubStep( const IndexType SubStepIndex, const LocalSystemVectorType& rFixedDofsValues, LocalSystemMatrixType& rIntermediateStepResidualVectors)`

# Calls

- [GetDeltaTime](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime.md)
- [GetIntegrationTheta](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetIntegrationTheta.md)
- [GetMatrixRowBegin](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowBegin.md)
- [GetMatrixRowEnd](../../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRowEnd.md)

# Called by

- [SolveWithLumpedMassMatrix](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix.md)