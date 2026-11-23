---
type: C++ Method
title: GetActualizeRHSflag
resource: solving_strategies/convergencecriterias/convergence_criteria.h#L227-L230
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/convergencecriterias/and_criteria/And_Criteria/And_Criteria
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/or_criteria/Or_Criteria/Or_Criteria
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/CyclicSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/IsConverged
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/IsConverged
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GetActualizeRHSflag()`

# Called by

- [And_Criteria](../../../../../functions/solving_strategies/convergencecriterias/and_criteria/And_Criteria/And_Criteria.md)
- [Or_Criteria](../../../../../functions/solving_strategies/convergencecriterias/or_criteria/Or_Criteria/Or_Criteria.md)
- [CyclicSolve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/CyclicSolve.md)
- [Solve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve.md)
- [IsConverged](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/IsConverged.md)
- [SolveSolutionStep](../../../../../functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/SolveSolutionStep.md)
- [IsConverged](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/IsConverged.md)
- [InitializeSolutionStep](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/InitializeSolutionStep.md)
- [SolveSolutionStep](../../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/SolveSolutionStep.md)