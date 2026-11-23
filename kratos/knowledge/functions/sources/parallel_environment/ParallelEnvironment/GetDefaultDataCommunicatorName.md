---
type: C++ Method
title: GetDefaultDataCommunicatorName
resource: sources/parallel_environment.cpp#L135-L139
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_data_communicator_factory/TestDataCommunicatorFactory/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_parallel_environment/TestParallelEnvironment/testDefaultRankAndSize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string ParallelEnvironment::GetDefaultDataCommunicatorName()`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE.md)
- [setUp](../../../../functions/mpi/tests/test_data_communicator_factory/TestDataCommunicatorFactory/setUp.md)
- [testDefaultRankAndSize](../../../../functions/mpi/tests/test_parallel_environment/TestParallelEnvironment/testDefaultRankAndSize.md)