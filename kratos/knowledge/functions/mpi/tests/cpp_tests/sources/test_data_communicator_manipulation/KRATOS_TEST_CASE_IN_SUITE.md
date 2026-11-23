---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/sources/test_data_communicator_manipulation.cpp#L29-L39
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/GetDefaultDataCommunicatorName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(CreateMPIDataCommunicatorDuplicate, KratosMPICoreFastSuite)`

# Calls

- [UnregisterDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicator.md)
- [HasDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicator.md)
- [GetDefaultDataCommunicatorName](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/GetDefaultDataCommunicatorName.md)
- [SetDefaultDataCommunicator](../../../../../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicator.md)