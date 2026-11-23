---
type: Python Method
title: ExecuteInitializeSolutionStep
resource: python_scripts/fix_scalar_variable_process.py#L48-L54
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/interval_utility/IntervalUtility/IsInInterval
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteInitializeSolutionStep(self):`

# Calls

- [IsInInterval](../../../../functions/utilities/interval_utility/IntervalUtility/IsInInterval.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)