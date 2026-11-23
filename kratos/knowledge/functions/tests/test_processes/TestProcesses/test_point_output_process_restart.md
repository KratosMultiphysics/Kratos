---
type: Python Method
title: test_point_output_process_restart
resource: tests/test_processes.py#L1542-L1608
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/SolutionLoopPointOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_point_output_process_restart(self):`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SolutionLoopPointOutputProcesses](../../../../functions/tests/test_processes/SolutionLoopPointOutputProcesses.md)