---
type: Python Method
title: _SetInitialTime
resource: tests/test_initial_time.py#L62-L67
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetVector
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_initial_time/TestInitialTime/testInitialTime
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _SetInitialTime(parameters, initial_time):`

# Calls

- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [SetVector](../../../../functions/sources/kratos_parameters/Parameters/SetVector.md)

# Called by

- [testInitialTime](../../../../functions/tests/test_initial_time/TestInitialTime/testInitialTime.md)