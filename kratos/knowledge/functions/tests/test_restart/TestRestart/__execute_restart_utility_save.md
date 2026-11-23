---
type: Python Method
title: __execute_restart_utility_save
resource: tests/test_restart.py#L172-L199
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/restart_utility/RestartUtility/IsRestartOutputStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/SaveRestart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_restart/TestRestart/test_restart_utility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __execute_restart_utility_save(self, model_part_name, restart_time): #creating a Model which will be destroyed at the end of the save function`

# Calls

- [IsRestartOutputStep](../../../../functions/python_scripts/restart_utility/RestartUtility/IsRestartOutputStep.md)
- [SaveRestart](../../../../functions/python_scripts/restart_utility/RestartUtility/SaveRestart.md)

# Called by

- [test_restart_utility](../../../../functions/tests/test_restart/TestRestart/test_restart_utility.md)