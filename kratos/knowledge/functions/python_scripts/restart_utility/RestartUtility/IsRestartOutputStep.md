---
type: Python Method
title: IsRestartOutputStep
resource: python_scripts/restart_utility.py#L157-L164
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/save_restart_process/SaveRestartProcess/IsOutputStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/__execute_restart_utility_save
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def IsRestartOutputStep(self):`

# Called by

- [IsOutputStep](../../../../functions/python_scripts/save_restart_process/SaveRestartProcess/IsOutputStep.md)
- [__execute_restart_utility_save](../../../../functions/tests/test_restart/TestRestart/__execute_restart_utility_save.md)