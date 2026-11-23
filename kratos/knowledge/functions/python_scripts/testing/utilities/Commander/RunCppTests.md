---
type: Python Method
title: RunCppTests
resource: python_scripts/testing/utilities.py#L237-L263
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/TestToAppName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/_RunTest
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/testing/run_cpp_tests/main
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunCppTests(self, applications, timer, config):`

# Calls

- [GetKratosMultiphysicsPath](../../../../../functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath.md)
- [TestToAppName](../../../../../functions/python_scripts/testing/utilities/Commander/TestToAppName.md)
- [_RunTest](../../../../../functions/python_scripts/testing/utilities/Commander/_RunTest.md)

# Called by

- [main](../../../../../functions/python_scripts/testing/run_cpp_tests/main.md)