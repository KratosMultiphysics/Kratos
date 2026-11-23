---
type: C++ Method
title: CyclicSolve
resource: solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h#L440-L621
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/MoveMeshFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/GetActualizeRHSflag
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CyclicSolve()`

# Calls

- [MoveMeshFlag](../../../../../functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/MoveMeshFlag.md)
- [GetActualizeRHSflag](../../../../../functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/GetActualizeRHSflag.md)

# Called by

- [Solve](../../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/Solve.md)