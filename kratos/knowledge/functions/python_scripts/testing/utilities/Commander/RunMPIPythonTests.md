---
type: Python Method
title: RunMPIPythonTests
resource: python_scripts/testing/utilities.py#L265-L308
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetListOfAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/_RunTest
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/testing/run_python_mpi_tests/main
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunMPIPythonTests(self, applications, mpi_command, mpi_flags, num_processes_flag, num_processes, level, verbose, command, timer):`

# Calls

- [GetListOfAvailableApplications](../../../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)
- [GetKratosMultiphysicsPath](../../../../../functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath.md)
- [_RunTest](../../../../../functions/python_scripts/testing/utilities/Commander/_RunTest.md)

# Called by

- [main](../../../../../functions/python_scripts/testing/run_python_mpi_tests/main.md)