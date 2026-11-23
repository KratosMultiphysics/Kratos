---
type: Python Method
title: __execute_restart_test
resource: tests/test_restart.py#L164-L170
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_restart/TestRestart/__execute_restart_save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/__execute_restart_load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/_check_modelpart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_restart/TestRestart/test_restart_NOTRACE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/test_restart_TRACE_ERROR
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/test_restart_TRACE_ALL
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __execute_restart_test(self, serializer_flag):`

# Calls

- [__execute_restart_save](../../../../functions/tests/test_restart/TestRestart/__execute_restart_save.md)
- [__execute_restart_load](../../../../functions/tests/test_restart/TestRestart/__execute_restart_load.md)
- [_check_modelpart](../../../../functions/tests/test_restart/TestRestart/_check_modelpart.md)

# Called by

- [test_restart_NOTRACE](../../../../functions/tests/test_restart/TestRestart/test_restart_NOTRACE.md)
- [test_restart_TRACE_ERROR](../../../../functions/tests/test_restart/TestRestart/test_restart_TRACE_ERROR.md)
- [test_restart_TRACE_ALL](../../../../functions/tests/test_restart/TestRestart/test_restart_TRACE_ALL.md)