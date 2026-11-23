---
type: C++ Function
title: Initialize
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L672-L702
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/schemes/scheme/Scheme/SchemeIsInitialized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/scheme/Scheme/ElementsAreInitialized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/scheme/Scheme/InitializeElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/scheme/Scheme/ConditionsAreInitialized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/scheme/Scheme/InitializeConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Initialize() override`

# Calls

- [SchemeIsInitialized](../../../../functions/solving_strategies/schemes/scheme/Scheme/SchemeIsInitialized.md)
- [ElementsAreInitialized](../../../../functions/solving_strategies/schemes/scheme/Scheme/ElementsAreInitialized.md)
- [InitializeElements](../../../../functions/solving_strategies/schemes/scheme/Scheme/InitializeElements.md)
- [ConditionsAreInitialized](../../../../functions/solving_strategies/schemes/scheme/Scheme/ConditionsAreInitialized.md)
- [InitializeConditions](../../../../functions/solving_strategies/schemes/scheme/Scheme/InitializeConditions.md)