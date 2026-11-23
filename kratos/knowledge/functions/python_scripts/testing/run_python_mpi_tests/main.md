---
type: Python Function
title: main
resource: python_scripts/testing/run_python_mpi_tests.py#L15-L94
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
  - target: functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/PrintTestSummary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def main():`

# Calls

- [GetPython3Command](../../../../functions/python_scripts/testing/utilities/GetPython3Command.md)
- [GetListOfAvailableApplications](../../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)
- [RunMPIPythonTests](../../../../functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests.md)
- [PrintTestSummary](../../../../functions/python_scripts/testing/utilities/PrintTestSummary.md)