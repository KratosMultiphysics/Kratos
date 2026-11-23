---
type: C++ Method
title: SubstepCount
resource: solving_strategies/strategies/butcher_tableau.h#L83
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaLastSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static constexpr unsigned int SubstepCount()`

# Called by

- [SolveSubStep](../../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)
- [SolveWithLumpedMassMatrix](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/SolveWithLumpedMassMatrix.md)
- [PerformRungeKuttaLastSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaLastSubStep.md)