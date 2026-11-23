---
type: C++ Function
title: SolveSolutionStep
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L919-L1105
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetCurrentSolution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/MoveMeshFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/GetActualizeRHSflag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetKeepSystemConstantDuringIterations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SolveSolutionStep() override`

# Calls

- [GetCurrentSolution](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetCurrentSolution.md)
- [MoveMeshFlag](../../../../functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/MoveMeshFlag.md)
- [GetActualizeRHSflag](../../../../functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/GetActualizeRHSflag.md)
- [GetKeepSystemConstantDuringIterations](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetKeepSystemConstantDuringIterations.md)