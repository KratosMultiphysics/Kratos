---
type: Python Method
title: test_save_restart_process
resource: tests/test_restart.py#L245-L296
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_restart/IsRestartFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/__execute_restart_load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/_check_modelpart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_save_restart_process(self):`

# Calls

- [IsRestartFile](../../../../functions/tests/test_restart/IsRestartFile.md)
- [__execute_restart_load](../../../../functions/tests/test_restart/TestRestart/__execute_restart_load.md)
- [_check_modelpart](../../../../functions/tests/test_restart/TestRestart/_check_modelpart.md)