---
type: Python Method
title: __execute_restart_load
resource: tests/test_restart.py#L154-L162
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_restart/TestRestart/__execute_restart_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/test_save_restart_process
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __execute_restart_load(self, current_model, file_name, serializer_flag):`

# Called by

- [__execute_restart_test](../../../../functions/tests/test_restart/TestRestart/__execute_restart_test.md)
- [test_save_restart_process](../../../../functions/tests/test_restart/TestRestart/test_save_restart_process.md)