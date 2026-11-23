---
type: Python Function
title: GetListOfAvailableApplications
resource: python_scripts/kratos_utilities.py#L53-L63
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_interface/python_registry/register_wrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/CheckIfApplicationsAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GetNotAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_cpp_mpi_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_cpp_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_python_mpi_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_python_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunPythonTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetListOfAvailableApplications():`

# Calls

- [GetKratosMultiphysicsPath](../../../functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath.md)

# Called by

- [register_wrapper](../../../functions/python_interface/python_registry/register_wrapper.md)
- [CheckIfApplicationsAvailable](../../../functions/python_scripts/kratos_utilities/CheckIfApplicationsAvailable.md)
- [GetNotAvailableApplications](../../../functions/python_scripts/kratos_utilities/GetNotAvailableApplications.md)
- [main](../../../functions/python_scripts/testing/run_cpp_mpi_tests/main.md)
- [main](../../../functions/python_scripts/testing/run_cpp_tests/main.md)
- [main](../../../functions/python_scripts/testing/run_python_mpi_tests/main.md)
- [main](../../../functions/python_scripts/testing/run_python_tests/main.md)
- [RunPythonTests](../../../functions/python_scripts/testing/utilities/Commander/RunPythonTests.md)
- [RunMPIPythonTests](../../../functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests.md)