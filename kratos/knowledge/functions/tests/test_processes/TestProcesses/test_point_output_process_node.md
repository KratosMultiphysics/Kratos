---
type: Python Method
title: test_point_output_process_node
resource: tests/test_processes.py#L1311-L1391
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/SolutionLoopPointOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_point_output_process_node(self):`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SolutionLoopPointOutputProcesses](../../../../functions/tests/test_processes/SolutionLoopPointOutputProcesses.md)
- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)