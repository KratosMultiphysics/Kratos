---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/sources/test_mpi_communicator.cpp#L75-L104
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/ModelPartForMPICommunicatorTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quaternion/W-2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetEquationId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorCreation, KratosMPICoreFastSuite)`

# Calls

- [MyPID](../../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [TotalProcesses](../../../../../../functions/sources/communicator/Communicator/TotalProcesses.md)
- [ModelPartForMPICommunicatorTests](../../../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/ModelPartForMPICommunicatorTests.md)
- [IsDefined](../../../../../../functions/containers/flags/IsDefined.md)
- [W](../../../../../../functions/utilities/quaternion/W-2.md)
- [SetEquationId](../../../../../../functions/includes/dof/Dof/SetEquationId.md)
- [GetNumberOfColors](../../../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [GlobalNumberOfNodes](../../../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)