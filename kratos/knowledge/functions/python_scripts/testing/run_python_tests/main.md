---
type: Python Function
title: main
resource: python_scripts/testing/run_python_tests.py#L9-L69
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
  - target: functions/python_scripts/testing/utilities/Commander/RunPythonTests
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
- [RunPythonTests](../../../../functions/python_scripts/testing/utilities/Commander/RunPythonTests.md)
- [PrintTestSummary](../../../../functions/python_scripts/testing/utilities/PrintTestSummary.md)