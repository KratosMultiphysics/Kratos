---
type: C++ Method
title: Solve
resource: solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h#L623-L729
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/CyclicSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/ReductionStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/PrintIncreaseOfDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/GetActualizeRHSflag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/scheme/Scheme/Clean
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double Solve() override`

# Calls

- [CyclicSolve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/CyclicSolve.md)
- [ReductionStrategy](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/ReductionStrategy.md)
- [PrintIncreaseOfDeltaTime](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/PrintIncreaseOfDeltaTime.md)
- [GetActualizeRHSflag](../../../../../functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/GetActualizeRHSflag.md)
- [Clean](../../../../../functions/solving_strategies/schemes/scheme/Scheme/Clean.md)