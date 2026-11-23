---
type: C++ Method
title: InitializeConstraints
resource: solving_strategies/schemes/scheme.h#L260-L266
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/scheme/Scheme/SetConstraintsAreInitialized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void InitializeConstraints(ModelPart& rModelPart)`

# Calls

- [SetConstraintsAreInitialized](../../../../../functions/solving_strategies/schemes/scheme/Scheme/SetConstraintsAreInitialized.md)

# Called by

- [Initialize](../../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/Initialize.md)