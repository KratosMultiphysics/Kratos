---
type: Python Function
title: GetKratosMultiphysicsPath
resource: python_scripts/kratos_utilities.py#L48-L51
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/kratos_utilities/GetListOfAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/IsMPIAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_cpp_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunPythonTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunCppTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunMPICppTests
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetKratosMultiphysicsPath():`

# Called by

- [GetListOfAvailableApplications](../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)
- [IsMPIAvailable](../../../functions/python_scripts/kratos_utilities/IsMPIAvailable.md)
- [main](../../../functions/python_scripts/testing/run_cpp_tests/main.md)
- [RunPythonTests](../../../functions/python_scripts/testing/utilities/Commander/RunPythonTests.md)
- [RunCppTests](../../../functions/python_scripts/testing/utilities/Commander/RunCppTests.md)
- [RunMPIPythonTests](../../../functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests.md)
- [RunMPICppTests](../../../functions/python_scripts/testing/utilities/Commander/RunMPICppTests.md)