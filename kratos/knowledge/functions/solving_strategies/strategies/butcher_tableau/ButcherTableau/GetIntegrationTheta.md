---
type: C++ Method
title: GetIntegrationTheta
resource: solving_strategies/strategies/butcher_tableau.h#L139-L143
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
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

`constexpr double GetIntegrationTheta(const unsigned int SubStepIndex) const`

# Called by

- [SolveSubStep](../../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)
- [PerformRungeKuttaIntermediateSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaIntermediateSubStep.md)
- [PerformRungeKuttaLastSubStep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/PerformRungeKuttaLastSubStep.md)