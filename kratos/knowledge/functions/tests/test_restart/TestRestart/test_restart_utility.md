---
type: Python Method
title: test_restart_utility
resource: tests/test_restart.py#L231-L243
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_restart/TestRestart/__execute_restart_utility_save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/__execute_restart_utility_load
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/_check_modelpart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_restart_utility(self): # Here we only test SERIALIZER_NO_TRACE since the others are tested in the simple tests`

# Calls

- [__execute_restart_utility_save](../../../../functions/tests/test_restart/TestRestart/__execute_restart_utility_save.md)
- [__execute_restart_utility_load](../../../../functions/tests/test_restart/TestRestart/__execute_restart_utility_load.md)
- [_check_modelpart](../../../../functions/tests/test_restart/TestRestart/_check_modelpart.md)