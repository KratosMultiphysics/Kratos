---
type: Python Method
title: __execute_restart_utility_load
resource: tests/test_restart.py#L201-L218
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/LoadRestart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_restart/TestRestart/test_restart_utility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __execute_restart_utility_load(self, current_model, model_part_name, restart_time):`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [LoadRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/LoadRestart.md)

# Called by

- [test_restart_utility](../../../../functions/tests/test_restart/TestRestart/test_restart_utility.md)