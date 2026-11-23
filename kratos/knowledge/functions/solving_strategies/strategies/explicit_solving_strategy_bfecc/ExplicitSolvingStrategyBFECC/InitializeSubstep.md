---
type: C++ Method
title: InitializeSubstep
resource: solving_strategies/strategies/explicit_solving_strategy_bfecc.h#L533-L561
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/CopySolutionStepData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/PerformSubstep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SubstepData InitializeSubstep(const Substep SubstepType)`

# Calls

- [CopySolutionStepData](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/CopySolutionStepData.md)

# Called by

- [PerformSubstep](../../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/PerformSubstep.md)