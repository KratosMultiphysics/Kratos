---
type: C++ Method
title: UnregisterDataCommunicator
resource: sources/parallel_environment.cpp#L123-L127
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicatorDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/utilities/test_pointer_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_data_communicator_factory/TestDataCommunicatorFactory/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/TearDown
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::UnregisterDataCommunicator(const std::string& Name)`

# Calls

- [UnregisterDataCommunicatorDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/UnregisterDataCommunicatorDetail.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_pointer_communicator/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE.md)
- [tearDown](../../../../functions/mpi/tests/test_data_communicator_factory/TestDataCommunicatorFactory/tearDown.md)
- [TearDown](../../../../functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/TearDown.md)