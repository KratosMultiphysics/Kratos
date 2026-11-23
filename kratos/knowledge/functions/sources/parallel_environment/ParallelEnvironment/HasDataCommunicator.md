---
type: C++ Method
title: HasDataCommunicator
resource: sources/parallel_environment.cpp#L129-L133
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicatorDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_environment/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_parallel_environment/TestParallelEnvironment/testDataCommunicatorRetrievalFromParallelEnvironment
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/SetUp
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ParallelEnvironment::HasDataCommunicator(const std::string& rName)`

# Calls

- [HasDataCommunicatorDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/HasDataCommunicatorDetail.md)

# Called by

- [InitializeMPIParallelRun](../../../../functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_data_communicator_manipulation/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_environment/KRATOS_TEST_CASE_IN_SUITE.md)
- [testDataCommunicatorRetrievalFromParallelEnvironment](../../../../functions/mpi/tests/test_parallel_environment/TestParallelEnvironment/testDataCommunicatorRetrievalFromParallelEnvironment.md)
- [SetUp](../../../../functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/SetUp.md)