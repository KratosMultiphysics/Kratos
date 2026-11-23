---
type: C++ Method
title: SetDefaultDataCommunicator
resource: sources/parallel_environment.cpp#L38-L42
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicatorDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_data_communicator_factory/TestDataCommunicatorFactory/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_parallel_environment/TestParallelEnvironment/testDefaultRankAndSize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::SetDefaultDataCommunicator(const std::string& rName)`

# Calls

- [SetDefaultDataCommunicatorDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetDefaultDataCommunicatorDetail.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE.md)
- [tearDown](../../../../functions/mpi/tests/test_data_communicator_factory/TestDataCommunicatorFactory/tearDown.md)
- [testDefaultRankAndSize](../../../../functions/mpi/tests/test_parallel_environment/TestParallelEnvironment/testDefaultRankAndSize.md)