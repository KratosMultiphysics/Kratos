---
type: Python Method
title: _SetModelPartName
resource: tests/test_initial_time.py#L44-L59
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_initial_time/TestInitialTime/testInitialTime
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _SetModelPartName(parameters, model_part_name):`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [testInitialTime](../../../../functions/tests/test_initial_time/TestInitialTime/testInitialTime.md)