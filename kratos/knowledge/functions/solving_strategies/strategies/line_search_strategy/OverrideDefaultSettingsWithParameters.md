---
type: C++ Function
title: OverrideDefaultSettingsWithParameters
resource: solving_strategies/strategies/line_search_strategy.h#L528-L538
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/line_search_strategy/LineSearchStrategy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void OverrideDefaultSettingsWithParameters( Parameters& rDefaultSettings, const double MaxIterations, const bool ReformDofSetAtEachStep, const bool CalculateReactions )`

# Calls

- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)

# Called by

- [LineSearchStrategy](../../../../functions/solving_strategies/strategies/line_search_strategy/LineSearchStrategy.md)