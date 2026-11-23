---
type: C++ Method
title: SolveWithLumpedMassMatrix
resource: solving_strategies/strategies/explicit_solving_strategy_bfecc.h#L309-L327
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/ExtractSolutionStepData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/PerformSubstep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/CorrectErrorAfterForwardsAndBackwards
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SolveWithLumpedMassMatrix() override`

# Calls

- [GetDeltaTime](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime.md)
- [ExtractSolutionStepData](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/ExtractSolutionStepData.md)
- [PerformSubstep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/PerformSubstep.md)
- [CorrectErrorAfterForwardsAndBackwards](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/CorrectErrorAfterForwardsAndBackwards.md)