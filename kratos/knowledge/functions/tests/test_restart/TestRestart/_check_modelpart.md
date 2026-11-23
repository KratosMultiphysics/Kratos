---
type: Python Method
title: _check_modelpart
resource: tests/test_restart.py#L38-L142
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfTables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_restart/TestRestart/__execute_restart_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/test_restart_utility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/test_save_restart_process
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _check_modelpart(self, model_part):`

# Calls

- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [NumberOfTables](../../../../functions/includes/model_part/NumberOfTables.md)
- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [__execute_restart_test](../../../../functions/tests/test_restart/TestRestart/__execute_restart_test.md)
- [test_restart_utility](../../../../functions/tests/test_restart/TestRestart/test_restart_utility.md)
- [test_save_restart_process](../../../../functions/tests/test_restart/TestRestart/test_save_restart_process.md)