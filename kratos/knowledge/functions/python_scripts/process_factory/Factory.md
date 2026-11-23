---
type: Python Function
title: Factory
resource: python_scripts/process_factory.py#L80-L91
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/TimeAveragingProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Factory(settings, Model):`

# Calls

- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [TimeAveragingProcess](../../../functions/processes/time_averaging_process/TimeAveragingProcess/TimeAveragingProcess.md)