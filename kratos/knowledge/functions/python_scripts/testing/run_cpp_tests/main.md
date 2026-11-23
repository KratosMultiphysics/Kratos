---
type: Python Function
title: main
resource: python_scripts/testing/run_cpp_tests.py#L10-L79
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/utilities/GetPython3Command
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GetListOfAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/remove
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunCppTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/PrintTestSummary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def main(): # Define the command`

# Calls

- [GetPython3Command](../../../../functions/python_scripts/testing/utilities/GetPython3Command.md)
- [GetListOfAvailableApplications](../../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)
- [remove](../../../../functions/sources/kratos_filesystem/remove.md)
- [GetKratosMultiphysicsPath](../../../../functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath.md)
- [RunCppTests](../../../../functions/python_scripts/testing/utilities/Commander/RunCppTests.md)
- [PrintTestSummary](../../../../functions/python_scripts/testing/utilities/PrintTestSummary.md)