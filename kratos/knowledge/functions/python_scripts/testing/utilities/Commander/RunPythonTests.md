---
type: Python Method
title: RunPythonTests
resource: python_scripts/testing/utilities.py#L172-L235
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
  - target: functions/python_scripts/testing/run_python_tests/main
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunPythonTests(self, applications, level, verbose, command, timer):`

# Calls

- [GetListOfAvailableApplications](../../../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)
- [GetKratosMultiphysicsPath](../../../../../functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath.md)
- [_RunTest](../../../../../functions/python_scripts/testing/utilities/Commander/_RunTest.md)

# Called by

- [main](../../../../../functions/python_scripts/testing/run_python_tests/main.md)