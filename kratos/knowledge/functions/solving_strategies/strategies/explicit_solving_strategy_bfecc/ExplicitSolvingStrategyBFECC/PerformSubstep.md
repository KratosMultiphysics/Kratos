---
type: C++ Method
title: PerformSubstep
resource: solving_strategies/strategies/explicit_solving_strategy_bfecc.h#L420-L471
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/InitializeSubstep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/SubstepData/TimeDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/FinalizeSubstep
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/SolveWithLumpedMassMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void PerformSubstep(const Substep SubstepType)`

# Calls

- [GetDeltaTime](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDeltaTime.md)
- [InitializeSubstep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/InitializeSubstep.md)
- [TimeDirection](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/SubstepData/TimeDirection.md)
- [FinalizeSubstep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/FinalizeSubstep.md)

# Called by

- [SolveWithLumpedMassMatrix](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/SolveWithLumpedMassMatrix.md)